#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
int parse_args(int, char **);

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


int main (int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_proc_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_procs_count);

  
  int ret = parse_args(argc, argv);
 
  /*
  if (ret)
  {
    return -1;
  }
  */

  /* Decompose points among processes */
  int p = num_points / world_procs_count;
  int q = num_points % world_procs_count;
  int proc_points_count = world_proc_rank < q ? p + 1 : p;
  int proc_point_start_idx = world_proc_rank * p + (world_proc_rank < q ? world_proc_rank : q);

  /* Decompose points among threads */
  p = proc_points_count / num_threads;
  q = proc_points_count % num_threads;
  int thread_points_counts[num_threads];
  int thread_point_start_idx[num_threads];
  int i;
  for (i = 0; i < num_threads; ++i)
  {
    thread_points_counts[i] = i < q ? p+1 : p;
    thread_point_start_idx[i] = (i * p + (i < q ? i : p));
  }

  /* Read points and centers from files */
  double *points = malloc(sizeof(double)*num_points*dim);
  double *centers = malloc(sizeof(double)*num_centers*dim);

  FILE *f = fopen(points_file, "rb");
  fread(points, sizeof(double), num_points*dim, f);
  fclose(f);

  f = fopen(centers_file, "rb");
  fread(centers, sizeof(double), num_centers*dim, f);
  fclose(f);

  /* Data structures for computation */
  int length_sums_and_counts = num_threads*num_centers*(dim+1);
  double thread_center_sums_and_counts[length_sums_and_counts];
  double proc_cluster_assignments[proc_points_count];

  int itr_count = 0;
  int converged = 0;

  /* Main computation loop */
  while (!converged && itr_count < max_iterations)
  {
    ++itr_count;
    reset_array(thread_center_sums_and_counts, length_sums_and_counts);

  }

  MPI_Finalize();

}

void reset_array(double *array, int length)
{
}

int parse_args(int argc, char **argv)
{
  int index;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "n:d:k:m:t:T:o::c:p:b:")) != -1)
    switch (c)
      {
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
        if (optopt == 'n' || optopt == 'd' ||  optopt == 'k' ||  optopt == 'm' ||  optopt == 't' ||  optopt == 'T' ||  optopt == 'c' ||  optopt == 'p' ||  optopt == 'b' ||  optopt == 'o')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
  if (world_proc_rank == 0)
  {
    printf ("Program Arguments\n");
    printf (" n = %d\n d = %d\n k = %d\n m = %d\n t = %lf\n T = %d\n o = %s\n c = %s\n p = %s\n b = %d\n",num_points, dim, num_centers, max_iterations, err_threshold, num_threads, output_file, centers_file, points_file, bind_threads);
  }

  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);
  return 0;
}
