#include <stdio.h>
#include <fstream>
#include <stdio.h>

//#define WIN32 1

//D:\Genome\Cram\DRR346006_1.fastq D:\Genome\Cram\DRR346006_1.fqb 150

const int blockSize = 10000000;

int main(int argc, char** argv) {
	int nRead, n32, readSize, ik, lastOk, L;
	int i32, j32, i32p, j32p, iNew;
	int* vOfor, *vOback;
	int lineStart[5], lineEnd[5], mOne[32];

	char inputFile[1000], outputFile[1000], * vData;
	char* vData_loc, * vDDD;

	long long fileSize, fPos, countR, fDiff, fPosBefore;

	bool isDone;

	FILE* fi, * fo;

	if (argc != 4) {
		printf("Usage:\n");
		printf("1) Input file\n");
		printf("2) Output file\n");
		printf("3) Length of reads\n\n");
		return -1;
	}

	vData = (char*)malloc(sizeof(char) * (blockSize + 2));

	sprintf(inputFile, "%s", argv[1]);
	sprintf(outputFile, "%s", argv[2]);
	nRead = atoi(argv[3]);

	n32 = nRead / 32;
	if (nRead % 32 > 0) n32++;

	vOfor = (int*)malloc(sizeof(int) * 4 * n32);
	vOback = (int*)malloc(sizeof(int) * 4 * n32);

	printf("Input file: \"%s\"\n", inputFile);
	printf("Output file: \"%s\"\n", outputFile);
	printf("Length of reads: %i\n", nRead);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -2;
	}

#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);
#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	fseeko64(fi, 0, SEEK_SET);
#endif

	fo = fopen(outputFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -3;
	}

	countR = 0;

	fwrite(&nRead, sizeof(int), 1, fo);
	fwrite(&countR, sizeof(long long), 1, fo);

	fPos = 0;

	for (int i = 0; i < 32; i++) {
		mOne[i] = 1 << i;
	}

	fPosBefore = -1;

	do {
		if (fPosBefore == fPos) break;
		fPosBefore = fPos;
		fDiff = fileSize - fPos;
		if (fDiff < 10) break;

#ifdef WIN32
		_fseeki64(fi, fPos, SEEK_SET);
#else
		fseeko64(fi, fPos, SEEK_SET);
#endif

		if (fDiff > blockSize) {
			fDiff = blockSize;
			fread(vData, sizeof(char), fDiff, fi);
		}
		else {
			fread(vData, sizeof(char), fDiff, fi);
			vData[fDiff] = '\n';
			vData[fDiff + 1] = ' ';
			fDiff += 2;
		}
		
		readSize = fDiff;
		vData_loc = vData;
		while (vData_loc[0] == '\n' || vData_loc[0] == '\r') {
			vData_loc++;
			readSize--;
		}
		fPos += (long long)(fDiff - readSize);

		lastOk = 0;

		for (;;) {
			lineStart[0] = lastOk;
			
			for (int k = 0; k < 4; k++) {
				isDone = false;
				for (int i = lineStart[k]; i < readSize; i++) {
					if (vData_loc[i] == '\n' || vData_loc[i] == '\r') {
						lineEnd[k] = i;
						isDone = true;
						break;
					}
				}
				if (!isDone) break;
				isDone = false;
				ik = lineEnd[k];
				for (int i = ik; i < readSize; i++) {
					if (!(vData_loc[i] == '\n' || vData_loc[i] == '\r')) {
						lineStart[k + 1] = i;
						isDone = true;
						break;
					}
				}
				
				if (!isDone)break;
			}

			if (!isDone) break;

			fPos += (lineStart[4] - lineStart[0]);

			lastOk = lineStart[4];

			for (int i = 0; i < 4 * n32; i++) {
				vOfor[i] = 0;
				vOback[i] = 0;
			}

			L = lineEnd[1] - lineStart[1];
			if (L >= 32 * n32) L = 32 * n32;
			vDDD = vData_loc + lineStart[1];

			for (int ii = 0; ii < L; ii++) {
				iNew = L - ii - 1;
				i32 = ii / 32;
				j32 = ii % 32;
				i32p = iNew / 32;
				j32p = iNew % 32;
				switch (vDDD[ii]) {
				case 'A':
				case 'a':
					vOfor[4 * i32 + 0] |= mOne[j32];
					vOback[4 * i32p + 3] |= mOne[j32p];
					break;
				case 'C':
				case 'c':
					vOfor[4 * i32 + 1] |= mOne[j32];
					vOback[4 * i32p + 2] |= mOne[j32p];
					break;
				case 'G':
				case 'g':
					vOfor[4 * i32 + 2] |= mOne[j32];
					vOback[4 * i32p + 1] |= mOne[j32p];
					break;
				case 'T':
				case 't':
					vOfor[4 * i32 + 3] |= mOne[j32];
					vOback[4 * i32p + 0] |= mOne[j32p];
					break;
				}
			}
			fwrite(vOfor, sizeof(int), 4 * n32, fo);
			fwrite(vOback, sizeof(int), 4 * n32, fo);
			countR++;
		}

	} while (true);

	free(vData); vData = nullptr;
	free(vOfor); vOfor = nullptr;
	free(vOback); vOback = nullptr;

#ifdef WIN32
	_fseeki64(fo, 4, SEEK_SET);
#else
	fseeko64(fo, 4, SEEK_SET);
#endif

	fwrite(&countR, sizeof(long long), 1, fo);

	fclose(fo); fo = nullptr;

	fclose(fi); fi = nullptr;
	return 0;
}
