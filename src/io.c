#include<stdio.h>
#include<stdint.h>

void write_to_file();

int main2(int argc, char *argv[])
{
//	char *file = "E:\\Sali\\git\\github\\esaliya\\ccpp\\KMeansC\\data\\100n_10k\\points.bin";
	//char *file = "E:\\Sali\\git\\github\\esaliya\\ccpp\\KMeansC\\data\\1n_1k\\points.bin";


	//write_to_file();
	char *file="/N/u/sekanaya/sali/git/github/esaliya/ccpp/KMeansC/data/1e6n_1000k/centers_LittleEndian.bin";
	FILE *f;
	f = fopen(file, "rb");
	char buffer[8];

	fread(&buffer, 8, 1, f);

	double x = *((double*)buffer);

	double y;
	fread(&y, 8, 1, f);
	
	printf("%0.9f\n", x);
	printf("%0.9f\n", y);
	
	int i;
	for (i = 0; i < 8; ++i)
	{
		printf("%x\n", buffer[i]);
	}

	fclose(f);

	getchar();
	return 0;
}

void write_to_file()
{
	char *file="/N/u/sekanaya/sali/git/github/esaliya/ccpp/KMeansC/data/1n_1k/points.bin";
	FILE *f;
	f = fopen(file, "wb");
	double x = 17;
	fwrite(&x, 8, 1, f);
	fclose(f);
}
