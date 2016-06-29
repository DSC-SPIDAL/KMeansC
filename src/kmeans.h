/*
 * kmeans.h
 *
 *  Created on: Jun 27, 2016
 *      Author: saliya
 */

#ifndef SRC_KMEANS_H_
#define SRC_KMEANS_H_

#include <time.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <sched.h>
#include <sys/syscall.h>


#endif /* SRC_KMEANS_H_ */

int parse_args(int argc, char **argv);
void set_thread_affinity(int world_proc_rank, int thread_id, int num_threads);
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
void get_lengths_array(int num_points, int procs_count, int *lengths);

/* MPI related variables */
char *machine_name;
int node_count;

int world_proc_rank;
int world_procs_count;

void print(const char* fmt, ...){
	if (world_proc_rank == 0){
		va_list args;
		va_start(args, fmt);
		vprintf(fmt, args);
		va_end(args);
	}
}

char* get_current_time(){
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	char* t = asctime (timeinfo);
	t[strlen(t)-1] = '\0';
	return t;
}

void print_affinity(int world_proc_rank, int thread_id){
	/* Print affinity */

	cpu_set_t mask;
	// We need the thread pid (even if we are in openmp)
	pid_t tid = (pid_t) syscall(SYS_gettid);
	// Get the affinity
	CPU_ZERO (&mask); // clear mask
	if (sched_getaffinity(tid, sizeof(mask), &mask) == -1) {
		printf("Error cannot do sched_getaffinity at rank %d and thread %d\n",
				world_proc_rank, thread_id);
	}

	char *bp;
	size_t size;
	FILE *stream;

	stream = open_memstream(&bp, &size);
	fprintf(stream, "Rank %d Thread %d, tid %d, affinity ", world_proc_rank,
			thread_id, tid);
	fflush(stream);

	// Print it
	int j;
	for (j = 0; j < CPU_SETSIZE; ++j) {
		if (CPU_ISSET(j, &mask)) {
			/*printf("Rank %d Thread %d, tid %d, affinity %d\n",
			 world_proc_rank, thread_id, tid, j);*/
			fprintf(stream, "%d ", j);
		}
	}
	fclose(stream);
	printf("%s\n", bp);
}
