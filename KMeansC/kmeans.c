#include<mpi.h>
#include<omp.h>
#include<stdio.h>
#include<stdint.h>
int main(int argc, char *argv[])
{
//	char *file = "E:\\Sali\\git\\github\\esaliya\\ccpp\\KMeansC\\data\\100n_10k\\points.bin";
	char *file = "E:\\Sali\\git\\github\\esaliya\\ccpp\\KMeansC\\data\\1n_1k\\points.bin";

	FILE *f;
	fopen_s(&f, file, "rb");
	char buffer[8];

	fread(&buffer, 8, 1, f);
	int i;
	for (i = 0; i < 8; ++i)
	{
		printf("%x\n", buffer[i]);
	}

	fclose(f);

	getchar();
	return 0;
}