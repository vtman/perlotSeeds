#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <chrono>

#include <omp.h>

//#define WIN32 1

const int NSIZE = 1000000;//Modify

//D:\Genome\library

class fileData {
public:
	fileData();
	~fileData();

	int allocateIndex();
	int processFile(int indFile);

	char inputFile[1000], outputFile[1000], inputTemplate[1000], outputTemplate[1000];
	int nn[16];
	char** vv, * v1;
	int* vIndex, nIndex;
	long long* vCount;
	int recordSize, nBits;
	int id, indexLevel;
	long long fileSize;
	FILE* fi, * fo;
	int nRecords;
};

int fileData::allocateIndex() {
	nIndex = 1 << indexLevel;
	vIndex = (int*)malloc(sizeof(int) * nIndex);
	vCount = (long long*)malloc(sizeof(long long) * NSIZE);
	for (int i = 0; i < NSIZE; i++) {
		vCount[i] = 0L;
	}
	return 0;
}

fileData::fileData() {
	vIndex = nullptr;
	vv = nullptr;
	v1 = nullptr;
	fi = nullptr;
	fo = nullptr;
}

fileData::~fileData() {
	if (vIndex != nullptr) { free(vIndex); vIndex = nullptr; }
	if (vCount != nullptr) { free(vCount); vCount = nullptr; }

	if (fi != nullptr) { fclose(fi); fi = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr; }
}

int fileData::processFile(int indFile) {
	int ii, jj, ival, uMask, start8, rem8;
	int indCurrent, nCount;
	int ntot;
	//printf("Thread: %i, file: %i\n", id, indFile);

	sprintf(inputFile, inputTemplate, indFile);
	sprintf(outputFile, outputTemplate, indFile);
	//printf("Input file: %s\n", inputFile);
	//printf("Output file: %s\n", outputFile);

	fi = fopen(inputFile, "rb");

#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
#else
	fseeko(fi, 0, SEEK_END);
	fileSize = ftello(fi);
#endif

	//printf("File size: %lli\n", fileSize);

	rewind(fi);
	nRecords = (int)(fileSize / ((long long)recordSize));
	//printf("Number of records: %i\n", nRecords);

	//printf("Thread: %i, file: %i, nRecords: %i\n", id, indFile, nRecords);

	v1 = (char*)malloc(sizeof(char) * nRecords * recordSize);

	vv = (char**)malloc(sizeof(char*) * 16);
	for (int i = 0; i < 16; i++) {
		vv[i] = (char*)malloc(sizeof(char) * nRecords * recordSize);
	}

	fread(v1, sizeof(char), nRecords * recordSize, fi);
	fclose(fi); fi = nullptr;


	uMask = 0;
	for (int i = 0; i < indexLevel; i++) {
		uMask |= (1 << i);
	}

	start8 = nBits - indexLevel;
	rem8 = start8 % 8;
	//if (rem8 > 0) start8 += 8;
	start8 /= 8;

	for (int i = 0; i < nIndex; i++) {
		vIndex[i] = 0;
	}

	for (int i = 0; i < nRecords; i++) {
		ival = uMask & (*(int*)(v1 + i * recordSize + start8) >> rem8);
		vIndex[ival]++;
	}
	for (int i = 1; i < nIndex; i++) {
		vIndex[i] += vIndex[i - 1];
	}

	for (int k = 0; k < nBits; k += 4) {
		ii = k / 8;
		jj = k % 8;

		for (int i = 0; i < 16; i++) {
			nn[i] = 0;
		}
		for (int i = 0; i < nRecords; i++) {
			ival = (*(unsigned char*)(v1 + i * recordSize + ii) >> jj) & 15;
			memcpy(vv[ival] + recordSize * nn[ival], v1 + recordSize * i, sizeof(char) * recordSize);
			nn[ival]++;

		}
		ntot = 0;
		for (int i = 0; i < 16; i++) {
			memcpy(v1 + ntot * recordSize, vv[i], sizeof(char) * nn[i] * recordSize);
			ntot += nn[i];
		}
		//memcpy(v1, v2, sizeof(char) * n2 * recordSize);
		//memcpy(v1 + n2 * recordSize, v3, sizeof(char) * n3 * recordSize);
	}

	indCurrent = 0;
	nCount = 1;
	for (int i = 1; i < nRecords; i++) {
		if (memcmp(v1 + indCurrent * recordSize, v1 + i * recordSize, (recordSize - 4) * sizeof(char)) == 0) {
			nCount++;
		}
		else {
			vCount[nCount]++;
			indCurrent = i;
			nCount = 1;
		}
	}
	vCount[nCount]++;

	fo = fopen(outputFile, "wb");
	fwrite(vIndex, sizeof(int), nIndex, fo);
	fwrite(v1, sizeof(char), recordSize * nRecords, fo);
	fclose(fo); fo = nullptr;

	for (int i = 0; i < 16; i++) {
		free(vv[i]);
	}
	free(vv); vv = nullptr;
	if (v1 != nullptr) { free(v1); v1 = nullptr; }

	return 0;
}

int main(int argc, char** argv) {
	int ik;
	char ioFolder[1000], inputFile[1000], inputTemplate[1000], outputTemplate[1000], sSeed[1000], sTemp[100];
	fileData** pfD;
	int seedLength, doubleWeight, nLetters, indexLevel;
	int nthreads, sizeRecord, nFiles, nBitsTotal, bits2process;
	FILE* fi;

	auto start = std::chrono::steady_clock::now();

	if (argc != 2) {
		printf("Usage:\n");
		printf("1) Input/output folder\n\n");
		return -1;
	}

	sprintf(ioFolder, "%s", argv[1]);

	printf("Input/output folder: %s\n\n", ioFolder);

	sprintf(inputFile, "%s/info.txt", ioFolder);
	printf("Input info file: %s\n", inputFile);
	fi = fopen(inputFile, "r");
	if (fi == nullptr) {
		printf("Error: cannot open file %s\n", inputFile);
		return -2;
	}
	fscanf(fi, "%s\t%i\n", sTemp, &seedLength);
	fscanf(fi, "%s\t%i\n", sTemp, &doubleWeight);
	fscanf(fi, "%s\t%s\n", sTemp, sSeed);
	fscanf(fi, "%s\t%i\n", sTemp, &nLetters);
	fscanf(fi, "%s\t%i\n", sTemp, &indexLevel);

	fclose(fi); fi = nullptr;

	printf("Seed length: %i\n", seedLength);
	printf("Double weight: %i\n", doubleWeight);
	printf("Seed: %s\n", sSeed);
	printf("Number of letters: %i\n", nLetters);
	printf("Index level: %i\n", indexLevel);

	bits2process = doubleWeight - 4 * nLetters;
	nBitsTotal = doubleWeight + 32 - 4 * nLetters;
	sizeRecord = nBitsTotal / 8;
	if (nBitsTotal % 8 > 0) sizeRecord++;
	nFiles = 1 << (4 * nLetters);

	//nFiles = 5000;

	printf("Bits to process: %i\n", bits2process);
	printf("Number of files: %i\n", nFiles);
	printf("Size of a record: %i bytes\n", sizeRecord);
	sprintf(inputTemplate, "%s/original/%%0%ix.bin", ioFolder, nLetters);
	sprintf(outputTemplate, "%s/sorted/%%0%ix.bin", ioFolder, nLetters);

	printf("Input template: %s\n", inputTemplate);
	printf("Output template: %s\n", outputTemplate);

#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}

	printf("Number of threads: %i\n", nthreads);

	pfD = (fileData**)malloc(sizeof(fileData*) * nthreads);

	for (int i = 0; i < nthreads; i++) {
		pfD[i] = new fileData();
		pfD[i]->id = i;
		sprintf(pfD[i]->inputTemplate, "%s", inputTemplate);
		sprintf(pfD[i]->outputTemplate, "%s", outputTemplate);
		pfD[i]->recordSize = sizeRecord;
		pfD[i]->nBits = bits2process;
		pfD[i]->indexLevel = indexLevel;
		pfD[i]->allocateIndex();
	}

	int ipercent, ipercentnew, nCount;
	ipercent = -1;
	nCount = 0;

#pragma omp parallel
	{
		int tid;
		fileData* fd;
		tid = omp_get_thread_num();
		fd = pfD[tid];
#pragma omp for
		for (ik = 0; ik < nFiles; ik++) {
			fd->processFile(ik);

#pragma omp critical
			{
				nCount++;
				ipercentnew = (int)(double(nCount) * 100.0 / double(nFiles));
				if (ipercentnew > ipercent) {
					ipercent = ipercentnew;
					printf("%i %%\n", ipercent);
				}
			}

		}

	}

	for (int k = 1; k < nthreads; k++) {
		for (int i = 0; i < NSIZE; i++) {
			pfD[0]->vCount[i] += pfD[k]->vCount[i];
		}
	}

	FILE* fStat;
	char outputFile[1000];
	sprintf(outputFile, "%s/stat.txt", ioFolder);
	fStat = fopen(outputFile, "w");
	for (int i = 0; i < NSIZE; i++) {
		if (pfD[0]->vCount[i] > 0) {
			fprintf(fStat, "%i\t%lli\n", i, pfD[0]->vCount[i]);
		}
	}
	fclose(fStat); fStat = nullptr;

	for (int i = 0; i < nthreads; i++) {
		if (pfD[i] != nullptr) {
			delete pfD[i]; pfD[i] = nullptr;
		}
	}
	free(pfD); pfD = nullptr;

	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
#ifdef WIN32
	printf("It took me %lld ms.\n\n", elapsed.count());
#else
	printf("It took me %ld ms.\n\n", elapsed.count());
#endif

	return 0;
}

