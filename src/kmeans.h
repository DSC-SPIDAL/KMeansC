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
