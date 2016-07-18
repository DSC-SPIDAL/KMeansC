#define _GNU_SOURCE

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <sched.h>
#include <sys/syscall.h>
#include "kmeans.h"
#include "juliet.h"

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

int verbose = 0;

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

	print("\nProgram Started on %s\n", get_current_time());

	print("  Reading points and centers ... ");
	double time = MPI_Wtime();
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
	print("\n    Done in %lf ms (on Rank 0)\n", (MPI_Wtime() - time)*1000);

	/* Data structures for computation */
	int length_sums_and_counts = num_threads * num_centers * (dim + 1);
	double *thread_centers_sums_and_counts = malloc(sizeof(double)*length_sums_and_counts);
	int *proc_clusters_assignments = malloc(sizeof(int)*proc_points_count);

	int itr_count = 0;
	int converged = 0;


	print("  Computing K-Means ... ");
	time = MPI_Wtime();
	/* Main computation loop */
	while (!converged && itr_count < max_iterations) {
		++itr_count;
		reset_array(thread_centers_sums_and_counts, length_sums_and_counts);

		if (num_threads > 1) {
			/* OpenMP parallel region */
			omp_set_num_threads(num_threads);
#pragma omp parallel
			{
				int effective_num_threads = omp_get_num_threads();
				int thread_id = omp_get_thread_num();
				if (effective_num_threads != num_threads && thread_id == 0){
					printf("Warning: Rank %d is running %d threads instead of the expected %d threads", world_proc_rank, effective_num_threads, num_threads);
				}

				if (bind_threads){
					set_thread_affinity(world_proc_rank, thread_id, num_threads);
				}

        
				if (verbose && itr_count == 1){
					print_affinity(world_proc_rank, thread_id);
				} 

				find_nearest_centers(points, centers, num_centers, dim,
									thread_centers_sums_and_counts, proc_clusters_assignments,
									thread_points_counts[thread_id], thread_points_start_idx[thread_id], thread_id*num_centers*(dim+1));

			}
		} else {
			if (bind_threads) {
				set_thread_affinity(world_proc_rank, 0, num_threads);
			}

			if (verbose && itr_count == 1) {
				print_affinity(world_proc_rank, 0);
			}

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
				// Note, can't break here as centers sums need to be divided to form new centers
				converged = 0;
			}

			for (d = 0; d < dim; ++d) {
				centers[(i * dim) + d] = thread_centers_sums_and_counts[(i
						* (dim + 1)) + d];
			}
		}
	}

  free(thread_centers_sums_and_counts);

	if (!converged) {
		if (world_proc_rank == 0) {
			printf(
					"\n      Stopping K-Means as max iteration count %d has reached",
					max_iterations);
		}
	}

	double times[1];
	times[0] = MPI_Wtime() - time;
	if (world_procs_count > 1){
		MPI_Reduce((world_proc_rank == 0 ? MPI_IN_PLACE : times), times, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}


	print("\n    Done in %d iterations and %lf ms (avg. across MPI)\n", itr_count, (times[0]*1000/world_procs_count));

	f = fopen(output_file, "w+");
	if (f != NULL) {
		int *recv = malloc(sizeof(int)*num_points);
        if (world_procs_count > 1) {
            // Gather cluster assignments
            print("  Gathering cluster assignments ...");
            time = MPI_Wtime();
            int lengths[world_procs_count];
            get_lengths_array(num_points, world_procs_count, lengths);
            int displas[world_procs_count];
            displas[0] = 0;
            for (i = 0; i < world_procs_count - 1; ++i){
            	displas[i+1] = lengths[i]+displas[i];
            }

            MPI_Allgatherv(proc_clusters_assignments, proc_points_count, MPI_INT, recv, lengths, displas, MPI_INT, MPI_COMM_WORLD);

            print("\n    Done in %lf ms (on Rank 0)\n", (MPI_Wtime() - time)*1000);
        }

        if (world_proc_rank == 0) {
        	double *all_points = malloc(sizeof(double) * num_points * dim);
        	FILE *fin = fopen(points_file, "rb");
        	fread(all_points, sizeof(double), num_points * dim, fin);
        	fclose(fin);

        	print("  Writing output file ...");
        	int *clusters = (world_procs_count > 1 ? recv : proc_clusters_assignments);

            time = MPI_Wtime();
            int d;
			for (i = 0; i < num_points; ++i) {
				fprintf(f, "%d\t", i);
				for (d = 0; d < dim; ++d){
					fprintf(f, "%lf\t", all_points[i*dim+d]);
				}
				fprintf(f, "%d\n", clusters[i]);
			}
            print("\n    Done in %lf ms (on Rank 0)\n", (MPI_Wtime() - time)*1000);
        }
        free(recv);
        free(proc_clusters_assignments);
        fclose(f);
	}

	free(points);
	free(centers);
	print("Program Terminated on %s\n", get_current_time());
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
		++centers_sums_and_counts[centers_offset + dim];
		accumulate(points, centers_sums_and_counts, points_offset,
				centers_offset, dim);
		clusters_assignments[i + points_start_idx] = min_dist_center;
	}

}

void accumulate(double *points, double *centers_sums_and_counts,
		int points_offset, int centers_offset, int dim) {
	int i;
	for (i = 0; i < dim; ++i) {
		centers_sums_and_counts[centers_offset + i] +=
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

void get_lengths_array(int num_points, int procs_count, int *lengths){
	int p = num_points / procs_count;
	int q = num_points % procs_count;
	int i;
	for (i = 0; i < procs_count; ++i) {
		lengths[i] = (i >= q ? p : p+1);
	}
}

void set_thread_affinity(int world_proc_rank, int thread_id, int num_threads){
	cpu_set_t mask;
	CPU_ZERO(&mask); // clear mask
	set_bit_mask(world_proc_rank, thread_id, num_threads, &mask);
	int ret = sched_setaffinity(0, sizeof(mask), &mask);
	if (ret < 0) {
		printf(
				"Error in setting thread affinity at rank %d and thread %d\n",
				world_proc_rank, thread_id);
	}
}

int parse_args(int argc, char **argv) {
	int index;
	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "n:d:k:m:t:T:o::c:p:b:v")) != -1)
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
		case 'v':
			verbose = 1;
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
