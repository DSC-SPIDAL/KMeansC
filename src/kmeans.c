#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <float.h>
#include <math.h>

int parse_args(int argc, char **argv);
void reset_array(double *array, int length);
double euclidean_distance(double *points1, double* points2, int offset1,
		int offset2, int dim);
int find_min_dist_center(double *points, double *centers, int num_centers,
		int dim, int points_offset);
void find_nearest_centers(double *points, double *centers, int num_centers,
		int dim, double *centers_sums_and_counts, int *clusters_assignments,
		int points_count, int points_start_idx, int offset);
void accumulate(double *points, double *centers_sums_and_counts,
		int points_offset, int centers_offset, int dim);

/* Configuration options */
int num_points;
int dim;
int num_centers;
int max_iterations;
double err_threshold;

char *points_file;
char *centers_file;
char *output_file;

int num_threads;
int bind_threads;

/* MPI related variables */
char *machine_name;
int node_count;

int world_proc_rank;
int world_procs_count;

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_proc_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_procs_count);

	int ret = parse_args(argc, argv);


	 if (ret) {
		 MPI_Finalize();
		 return -1;
	 }


	/* Decompose points among processes */
	int p = num_points / world_procs_count;
	int q = num_points % world_procs_count;
	int proc_points_count = world_proc_rank < q ? p + 1 : p;
	int proc_points_start_idx = world_proc_rank * p
			+ (world_proc_rank < q ? world_proc_rank : q);

	/* Decompose points among threads */
	p = proc_points_count / num_threads;
	q = proc_points_count % num_threads;
	int thread_points_counts[num_threads];
	int thread_points_start_idx[num_threads];
	int i;
	for (i = 0; i < num_threads; ++i) {
		thread_points_counts[i] = i < q ? p + 1 : p;
		thread_points_start_idx[i] = (i * p + (i < q ? i : q));
	}

	/* Read points and centers from files */
	double *points = malloc(sizeof(double) * proc_points_count * dim);
	double *centers = malloc(sizeof(double) * num_centers * dim);

	FILE *f = fopen(points_file, "rb");
	fseek(f, proc_points_start_idx * dim * sizeof(double), SEEK_SET);
	fread(points, sizeof(double), proc_points_count * dim, f);
	fclose(f);

	f = fopen(centers_file, "rb");
	fread(centers, sizeof(double), num_centers * dim, f);
	fclose(f);

	/* Data structures for computation */
	int length_sums_and_counts = num_threads * num_centers * (dim + 1);
	double thread_centers_sums_and_counts[length_sums_and_counts];
	int proc_clusters_assignments[proc_points_count];

	int itr_count = 0;
	int converged = 0;

	/* Main computation loop */
	while (!converged && itr_count < max_iterations) {
		++itr_count;
		reset_array(thread_centers_sums_and_counts, length_sums_and_counts);

		if (num_threads > 1) {
			/* OpenMP parallel region */
			// remember to send this offset threadIdx*numCenters*(dimension+1)
		} else {
			find_nearest_centers(points, centers, num_centers, dim,
					thread_centers_sums_and_counts, proc_clusters_assignments,
					thread_points_counts[0], thread_points_start_idx[0], 0);
		}

		if (num_threads > 1) {
			int t, c, d;
			for (t = 1; t < num_threads; ++t) {
				for (c = 0; c < num_centers; ++c) {
					for (d = 0; d < (dim + 1); ++d) {
						int offset_within_thread = (c * (dim + 1)) + d;
						thread_centers_sums_and_counts[offset_within_thread] +=
								thread_centers_sums_and_counts[(t * num_centers
										* (dim + 1)) + offset_within_thread];
					}
				}
			}
		}

		if (world_procs_count > 1) {
			MPI_Allreduce(MPI_IN_PLACE, thread_centers_sums_and_counts,
					num_centers * (dim + 1), MPI_DOUBLE, MPI_SUM,
					MPI_COMM_WORLD);
		}

		converged = 1;
		int d;
		double dist;
		for (i = 0; i < num_centers; ++i) {
			for (d = 0; d < dim; ++d) {
				thread_centers_sums_and_counts[(i * (dim + 1)) + d] /=
						thread_centers_sums_and_counts[(i * (dim + 1)) + dim];
			}

			dist = euclidean_distance(thread_centers_sums_and_counts, centers,
					(i * (dim + 1)), i * dim, dim);
			if (dist > err_threshold) {
				// Note, can't break heare as centers sums need to be divided to form new centers
				converged = 0;
			}

			for (d = 0; d < dim; ++d) {
				centers[(i * dim) + d] = thread_centers_sums_and_counts[(i
						* (dim + 1)) + d];
			}
		}
	}

	if (!converged) {
		if (world_proc_rank == 0) {
			printf(
					"    Stopping K-Means as max iteration count %d has reached\n",
					max_iterations);
		}
	}

	if (world_proc_rank == 0) {
		printf("    Done in %d iterations\n", itr_count);
	}

	MPI_Finalize();
	return 0;
}

void find_nearest_centers(double *points, double *centers, int num_centers,
		int dim, double *centers_sums_and_counts, int *clusters_assignments,
		int points_count, int points_start_idx, int offset) {
	int i;
	for (i = 0; i < points_count; ++i) {
		int points_offset = (points_start_idx + i) * dim;
		int min_dist_center = find_min_dist_center(points, centers, num_centers,
				dim, points_offset);
		int centers_offset = offset + min_dist_center * (dim + 1);
		++centers_sums_and_counts[centers_offset];
		accumulate(points, centers_sums_and_counts, points_offset,
				centers_offset, dim);
		clusters_assignments[i + points_start_idx] = min_dist_center;
	}

}

void accumulate(double *points, double *centers_sums_and_counts,
		int points_offset, int centers_offset, int dim) {
	int i;
	for (i = 0; i < dim; ++i) {
		centers_sums_and_counts[centers_offset + 1] +=
				points[points_offset + i];
	}
}

int find_min_dist_center(double *points, double *centers, int num_centers,
		int dim, int points_offset) {
	double min_dist = DBL_MAX;
	int min_dist_idx = -1;
	int i;
	for (i = 0; i < num_centers; ++i) {
		double dist = euclidean_distance(points, centers, points_offset,
				i * dim, dim);
		if (dist < min_dist) {
			min_dist = dist;
			min_dist_idx = i;
		}
	}
	return min_dist_idx;
}

double euclidean_distance(double *points1, double* points2, int offset1,
		int offset2, int dim) {
	double d = 0.0;
	double tmp;
	int i;
	for (i = 0; i < dim; ++i) {
		tmp = points1[i + offset1] - points2[i + offset2];
		d += tmp * tmp;
	}
	return sqrt(d);
}

void reset_array(double *array, int length) {
	int i;
	for (i = 0; i < length; ++i) {
		array[i] = 0.0;
	}
}

int parse_args(int argc, char **argv) {
	int index;
	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "n:d:k:m:t:T:o::c:p:b:")) != -1)
		switch (c) {
		case 'n':
			num_points = atoi(optarg);
			break;
		case 'd':
			dim = atoi(optarg);
			break;
		case 'k':
			num_centers = atoi(optarg);
			break;
		case 'm':
			max_iterations = atoi(optarg);
			break;
		case 't':
			err_threshold = atof(optarg);
			break;
		case 'T':
			num_threads = atoi(optarg);
			break;
		case 'o':
			output_file = optarg;
			break;
		case 'c':
			centers_file = optarg;
			break;
		case 'p':
			points_file = optarg;
			break;
		case 'b':
			bind_threads = atoi(optarg);
			break;
		case '?':
			if (optopt == 'n' || optopt == 'd' || optopt == 'k' || optopt == 'm'
					|| optopt == 't' || optopt == 'T' || optopt == 'c'
					|| optopt == 'p' || optopt == 'b' || optopt == 'o')
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint(optopt))
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			abort();
		}
	if (world_proc_rank == 0) {
		printf("Program Arguments\n");
		printf(
				" n = %d\n d = %d\n k = %d\n m = %d\n t = %lf\n T = %d\n o = %s\n c = %s\n p = %s\n b = %d\n",
				num_points, dim, num_centers, max_iterations, err_threshold,
				num_threads, output_file, centers_file, points_file,
				bind_threads);
	}

	for (index = optind; index < argc; index++)
		printf("Non-option argument %s\n", argv[index]);
	return 0;
}
