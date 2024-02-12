#include <fstream>
#include <stdlib.h>
#include <stdio.h>

//#define WIN32 1

//D:\Genome\Cram\GCA_000001405.29_GRCh38.p14_genomic.fna D:\Genome\Cram\ref38 200

//C:\Users\Sofya\source\repos\Gtools\Gtools\chr4.fna D:\Genome\Library4\chr4 200

const int N = 300000000;
const int M = 100000;

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], outputPrefix[1000];
	char mString[1000];
	FILE* fi, * foInfo, * foBin;
	char* vi;
	int* vo;
	int nGap, nCount, iGap32;
	int ipos1, ipos2, ipos3, nb, nLoc;
	int mOne[32], * iBlockGap;
	bool isFound, isRepeat;

	long long fileSize, fpos, chunkSize, oPosition, fposPrev;

	for (int i = 0; i < 32; i++) {
		mOne[i] = 1 << i;
	}

	oPosition = 0;

	fi = nullptr;
	foInfo = nullptr;
	foBin = nullptr;

	if (argc != 4) {
		printf("Error: wrong number of input parameters.\n");
		printf("1) path to FNA file\n");
		printf("2) output folder with prefix\n");
		printf("3) size of a gap\n\n");
		return -1;
	}

	sprintf(inputFile, "%s", argv[1]);
	sprintf(outputPrefix, "%s", argv[2]);
	nGap = atoi(argv[3]);
	printf("Input file: \"%s\".\n", inputFile);
	printf("Output path: \"%s\".\n", outputPrefix);
	printf("Gap: %i\n\n", nGap);

	iGap32 = nGap / 32;
	if (nGap % 32 > 0)iGap32++;
	iGap32 *= 4;

	iBlockGap = (int*)malloc(sizeof(int) * iGap32);
	for (int i = 0; i < iGap32; i++) {
		iBlockGap[i] = 0;
	}

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -2;
	}

	sprintf(outputFile, "%s.acgt", outputPrefix);
	foBin = fopen(outputFile, "wb");
	if (foBin == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -3;
	}

	sprintf(outputFile, "%s.iacgt", outputPrefix);
	foInfo = fopen(outputFile, "w");
	if (foInfo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -4;
	}

	vi = (char*)malloc(sizeof(char) * N);
	vo = (int*)malloc(4 * sizeof(int) * M);

#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);
#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
	fseeko64(fi, 0, SEEK_SET);
#endif

	printf("File size: %lld\n", fileSize);

	chunkSize = N;

	fposPrev = -1;

	do {
#ifdef WIN32
		fpos = _ftelli64(fi);
#else
		fpos = ftello64(fi);
#endif
		if (fpos == fposPrev) break;
		fposPrev = fpos;
		//printf("Posistion: %lld\n", fpos);
		if (fpos == fileSize) break;
		if (fpos + N > fileSize) {
			chunkSize = fileSize - fpos;
		}
		fread(vi, sizeof(char), chunkSize, fi);

		ipos1 = 0;
		do {
			isRepeat = false;
			ipos2 = ipos1 + 2;
			isFound = false;
			while (ipos2 < chunkSize) {
				if (vi[ipos2] == '>') {
					isFound = true;
					ipos2--;
					break;
				}
				ipos2++;
			}

			if (isFound) {
				if (vi[ipos1 + 1] == 'C' && vi[ipos1 + 2] == 'M') break;
				isRepeat = true;
				ipos1 = ipos2 + 1;
			}
			else {
				ipos2 = chunkSize - 1;
			}
		} while (isRepeat);

		if (vi[ipos1 + 1] != 'C' || vi[ipos1 + 2] != 'M') {
			fpos += (long long)ipos1;
#ifdef WIN32
			_fseeki64(fi, fpos, SEEK_SET);
#else
			fseeko64(fi, fpos, SEEK_SET);
#endif
			continue;
		}else {
			if (ipos1 != 0) {
				fpos += (long long)ipos1;
#ifdef WIN32
				_fseeki64(fi, fpos, SEEK_SET);
#else
				fseeko64(fi, fpos, SEEK_SET);
#endif
				continue;
			}
		}

		ipos3 = ipos1 + 2;
		while (vi[ipos3] != '\n' && vi[ipos3] != '\r') {
			ipos3++;
		}
		for (int i = 0; i < ipos3 - ipos1; i++) {
			mString[i] = vi[ipos1 + i];
		}
		mString[ipos3 - ipos1] = '\0';

		while (vi[ipos3] == '\n' && vi[ipos3] == '\r') {
			ipos3++;
		}

		fprintf(foInfo, "%s\n", mString);
		fprintf(foInfo, "%lli\n", oPosition);
		nCount = 0;
		nb = 0;
		nLoc = 0;
		vo[4 * nb + 0] = 0;
		vo[4 * nb + 1] = 0;
		vo[4 * nb + 2] = 0;
		vo[4 * nb + 3] = 0;
		for (int i = ipos3; i <= ipos2; i++) {
			if (vi[i] == '\n' || vi[i] == '\r') continue;
			switch (vi[i]) {
			case 'A':
			case 'a':
				vo[4 * nb + 0] |= mOne[nLoc];
				break;
			case 'C':
			case 'c':
				vo[4 * nb + 1] |= mOne[nLoc];
				break;
			case 'G':
			case 'g':
				vo[4 * nb + 2] |= mOne[nLoc];
				break;
			case 'T':
			case 't':
				vo[4 * nb + 3] |= mOne[nLoc];
				break;
			}
			nLoc++;
			if (nLoc == 32) {
				nb++;
				if (nb == M) {
					fwrite(vo, sizeof(int), 4 * nb, foBin);
					nb = 0;
				}
				vo[4 * nb + 0] = 0;
				vo[4 * nb + 1] = 0;
				vo[4 * nb + 2] = 0;
				vo[4 * nb + 3] = 0;

				nCount += 32;
				nLoc = 0;
			}
		}
		if (nLoc > 0) {
			nCount += nLoc;
			nb++;
		}
		if (nb > 0) {
			fwrite(vo, sizeof(int), 4 * nb, foBin);
		}
		
		fwrite(iBlockGap, sizeof(int), iGap32, foBin);
		
		fprintf(foInfo, "%i\n", nCount);

		oPosition += (iGap32 * 4);
		if ((nCount % 32) > 0) nCount += 32;
		nCount /= 32;
		oPosition += (nCount * 16);

		printf("XXX\t");
		for (int i = 0; i < 20; i++) {
			printf("%1c", vi[ipos1 + i]);
		}
		printf("\n");

		fpos += (long long)(ipos2 + 1);
#ifdef WIN32
		_fseeki64(fi, fpos, SEEK_SET);
#else
		fseeko64(fi, fpos, SEEK_SET);
#endif
	} while (true);

	free(vi); vi = nullptr;
	free(vo); vo = nullptr;

	fclose(fi); fi = nullptr;
	fclose(foInfo); foInfo = nullptr;
	fclose(foBin); foBin = nullptr;
	return 0;
}
