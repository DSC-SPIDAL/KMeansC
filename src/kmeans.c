#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int parse_args(int, char **);

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

int main (int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  
  int ret = parse_args(argc, argv);
 
  if (ret)
  {
    return -1;
  }


  MPI_Finalize();

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
  printf ("Program Arguments\n");
  printf (" n = %d\n d = %d\n k = %d\n m = %d\n t = %lf\n T = %d\n o = %s\n c = %s\n p = %s\n b = %d\n",num_points, dim, num_centers, max_iterations, err_threshold, num_threads, output_file, centers_file, points_file, bind_threads);

  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);
  return 0;
}
