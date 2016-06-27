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


#endif /* SRC_KMEANS_H_ */

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
