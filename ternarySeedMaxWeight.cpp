//E:\Temp2\TSeeds\Zen\TernaryBlocks\T3V4 E:\Temp2\perlotSeeds\ternarySeeds 3 4 3 50 50 80

#include <fstream>

const int MEM_BLOCK_SIZE = 10000000;

class findMW {
public:
	findMW();
	~findMW();

	FILE* fi, ** pfOut;

	char inputFolder[1000], outputFolder[1000], inputFile[1000], outputFile[1000];
	char sString[5000];
	char* vData;

	int mT, mV, blockSizeMin, blockSizeMax, readLenMin, readLenMax;
	int blockSM1, readLM1;
	int* wTT, * indTT;

	int* matCheck, * weightReadMax, * weightBlock, * countLen4Block, * indexLen4Block;

	long long fileSize;

	int findMatrix();
	int findBlockWeight();
	int allocateMemory();
	int printInfoParameters();
	int printInputValues();
	int readInputParameters(int na, char** va);
	int startProcessing(int na, char** va);
	int findMaxWeight();
	int reportMaxWeight();
};

findMW::findMW() {
	fi = nullptr;
	pfOut = nullptr;
	wTT = nullptr;
	indTT = nullptr;
	vData = nullptr;
	matCheck = nullptr;
	weightReadMax = nullptr;
	weightBlock = nullptr;
	countLen4Block = nullptr;
	indexLen4Block = nullptr;
}

findMW::~findMW() {
	if (fi != nullptr) { fclose(fi); fi = nullptr; }
	if (pfOut != nullptr) {
		for (int i = readLenMin; i <= readLenMax; i++) {
			if (pfOut[i] != nullptr) {
				fclose(pfOut[i]);
			}
		}
		free(pfOut); pfOut = nullptr;
	}
	if (countLen4Block != nullptr) { free(countLen4Block); countLen4Block = nullptr; }
	if (indexLen4Block != nullptr) { free(indexLen4Block); indexLen4Block = nullptr; }
	if (wTT != nullptr) { free(wTT); wTT = nullptr; }
	if (indTT != nullptr) { free(indTT); indTT = nullptr; }
	if (vData != nullptr) { free(vData); vData = nullptr; }
	if (matCheck != nullptr) { free(matCheck); matCheck = nullptr; }
	if (weightReadMax != nullptr) { free(weightReadMax); weightReadMax = nullptr; }
	if (weightBlock != nullptr) { free(weightBlock); weightBlock = nullptr; }
}

int findMW::startProcessing(int na, char** va) {
	if (na != 9) {
		printf("Error: wrong number of arguments\n");
		printInfoParameters();
		return -1;
	}
	if (readInputParameters(na, va) != 0) return -2;
	if (printInputValues() != 0) return -3;
	if (allocateMemory() != 0) return -4;
	if (findBlockWeight() != 0)return -5;
	if (findMatrix() != 0) return -6;
	if (findMaxWeight() != 0) return -7;
	if (findMatrix() != 0) return -8;
	if (reportMaxWeight() != 0) return -7;
	return 0;
}

int findMW::printInfoParameters() {
	printf("Input parameters:\n");
	printf("1) Input folder\n");
	printf("2) Output folder\n");
	printf("3) Number of transitions\n");
	printf("4) Number of transversions\n");
	printf("5) Block size (min)\n");
	printf("6) Block size (max)\n");
	printf("7) Read length (min)\n");
	printf("8) Read length (max)\n");
	return 0;
}

int findMW::readInputParameters(int na, char** va) {
	sprintf(inputFolder, "%s", va[1]);
	sprintf(outputFolder, "%s", va[2]);
	mT = atoi(va[3]);
	mV = atoi(va[4]);
	blockSizeMin = atoi(va[5]);
	blockSizeMax = atoi(va[6]);
	readLenMin = atoi(va[7]);
	readLenMax = atoi(va[8]);
	return 0;
}

int findMW::printInputValues() {
	printf("Input parameters:\n");
	printf("1) Input folder: %s\n", inputFolder);
	printf("2) Output folder: %s\n", outputFolder);
	printf("3) Number of transitions: %i\n", mT);
	printf("4) Number of transversions: %i\n", mV);
	printf("5) Block size (min): %i\n", blockSizeMin);
	printf("6) Block size (max): %i\n", blockSizeMax);
	printf("7) Read length (min): %i\n", readLenMin);
	printf("8) Read length (max): %i\n\n", readLenMax);
	return 0;
}

int findMW::allocateMemory() {
	readLM1 = readLenMax + 1;
	blockSM1 = blockSizeMax + 1;

	matCheck = (int*)malloc(sizeof(int) * readLM1 * blockSM1);
	for (int i = 0; i < readLM1 * blockSM1; i++) {
		matCheck[i] = 1;
	}

	weightReadMax = (int*)malloc(sizeof(int) * readLM1);
	for (int i = 0; i < readLM1; i++) {
		weightReadMax[i] = 0;
	}

	weightBlock = (int*)malloc(sizeof(int) * blockSM1);

	vData = (char*)malloc(sizeof(char) * MEM_BLOCK_SIZE);

	indTT = (int*)malloc(sizeof(int) * 2 * blockSM1);
	wTT = (int*)malloc(sizeof(int) * 2 * blockSM1);

	countLen4Block = (int*)malloc(sizeof(int) * blockSM1);
	indexLen4Block = (int*)malloc(sizeof(int) * blockSM1 * readLM1);
	return 0;
}

int findMW::findBlockWeight() {
	int rMin, n4, ms, sumTT, mB, mR;
	int bitAB, bitI;
	for (int ib = blockSizeMin; ib <= blockSizeMax; ib++) {
		weightBlock[ib] = 0;
		sprintf(inputFile, "%s/tblock_len_%i_mT_%i_mV_%i.bin", inputFolder, ib, mT, mV);
		printf("Input file (weight): \"%s\"\n", inputFile);
		fi = fopen(inputFile, "rb");

		if (fi == nullptr) {
			printf("Warning: cannot open the file\n");
			continue;
		}

		rMin = 2 * ib - 1;
		if (rMin < readLenMin) rMin = readLenMin;

		n4 = ib / 4;
		if (ib % 4 > 0) n4++;

		fread(vData, sizeof(char), n4, fi);

		fclose(fi); fi = nullptr;

		for (int j = 0; j < n4; j++) {
			bitI = vData[j];
			for (int k = 0; k < 4; k++) {
				bitAB = bitI & 3;
				indTT[4 * j + k] = (bitAB >> 1) + ((bitAB & 1) << 1);
				bitI = bitI >> 2;
			}
		}
		wTT[0] = 0;
		for (int j = 0; j < ib; j++) {
			wTT[j + 1] = wTT[j] + indTT[j];
		}
		for (int j = 0; j < ib; j++) {
			wTT[ib + j + 1] = wTT[ib + j] + indTT[j];
		}
		for (int j = 0; j < ib; j++) {
			indTT[j + ib] = indTT[j];
		}

		weightBlock[ib] = wTT[ib];

		for (int ir = rMin; ir <= readLenMax; ir++) {
			ms = ir - ib + 1;
			mB = ms / ib;
			mR = ms - mB * ib;
			if ((mB + 1) * wTT[ib] < weightReadMax[ir])continue;
			for (int ij = 0; ij < ib; ij++) {
				if (indTT[ij] == 0)continue;
				sumTT = mB * wTT[ib] + wTT[ij + mR] - wTT[ij];
				if (weightReadMax[ir] < sumTT)weightReadMax[ir] = sumTT;
			}
		}
	}

	for (int ib = blockSizeMin; ib <= blockSizeMax; ib++) {
		printf("Block weight for %i: %i\n", ib, weightBlock[ib]);
	}

	return 0;
}

int findMW::findMatrix() {
	int ms, mB, mR, val1, val2, nNonZeroElements;

	nNonZeroElements = 0;
	for (int ib = blockSizeMin; ib <= blockSizeMax; ib++) {
		if (weightBlock[ib] == 0) {
			for (int ir = readLenMin; ir <= readLenMax; ir++) {
				matCheck[ib * readLM1 + ir] = 0;
			}
			continue;
		}
		for (int ir = readLenMin; ir <= readLenMax; ir++) {
			if (matCheck[ib * readLM1 + ir] == 0)continue;
			ms = ir - ib + 1;
			if (ms <= 0) {
				matCheck[ib * readLM1 + ir] = 0;
				continue;
			}
			mB = ms / ib;
			if (mB < 1) {
				matCheck[ib * readLM1 + ir] = 0;
				continue;
			}
			mR = ms - mB * ib;
			val1 = weightBlock[ib] * mB + mR * 2;
			val2 = weightBlock[ib] * (mB + 1);
			if (val2 < val1) val1 = val2;
			if (weightReadMax[ir] > val1) {
				matCheck[ib * readLM1 + ir] = 0;
				continue;
			}
			nNonZeroElements++;
		}
	}

	for (int i = 0; i < blockSizeMin; i++) {
		countLen4Block[i] = 0;
	}
	for (int i = blockSizeMin; i <= blockSizeMax; i++) {
		countLen4Block[i] = 0;
		for (int ir = readLenMin; ir <= readLenMax; ir++) {
			if (matCheck[i * readLM1 + ir] == 0)continue;
			indexLen4Block[i * readLM1 + countLen4Block[i]] = ir;
			countLen4Block[i]++;
		}
	}

	printf("Number of non-zero elements of the matrix: %i\n", nNonZeroElements);
	return 0;
}

int findMW::findMaxWeight() {
	int n4, rowsPerBlock, bitAB, bitI, ir;
	int ms, mB, mR;
	long long rowsTotal;
	int nBlocks, nRem, nread, sumTT, sumQ;
	char* vLoc;
	for (int ib = blockSizeMin; ib <= blockSizeMax; ib++) {

		if (countLen4Block[ib] == 0) continue;

		sprintf(inputFile, "%s/tblock_len_%i_mT_%i_mV_%i.bin", inputFolder, ib, mT, mV);
		printf("Input file (first): \"%s\"\n", inputFile);
		fi = fopen(inputFile, "rb");

		if (fi == nullptr) {
			printf("Warning: cannot open the file\n");
			continue;
		}

		n4 = ib / 4;
		if (ib % 4 > 0) n4++;

		rowsPerBlock = MEM_BLOCK_SIZE / n4;

#ifdef WIN32
		_fseeki64(fi, 0, SEEK_END);
		fileSize = _ftelli64(fi);
#else
		fseeko64(fi, 0, SEEK_END);
		fileSize = ftello64(fi);
#endif

		rowsTotal = fileSize / ((long long)n4);
		printf("Rows: %lli\n", rowsTotal);

		nBlocks = (int)(rowsTotal / ((long long)rowsPerBlock));
		nRem = int(rowsTotal - ((long long)nBlocks) * ((long long)rowsPerBlock));

		if (nRem == 0) {
			nBlocks--;
			nRem = rowsPerBlock;
		}

#ifdef WIN32
		_fseeki64(fi, 0, SEEK_SET);
#else
		fseeko64(fi, 0, SEEK_SET);
#endif

		nread = rowsPerBlock;

		for (int nc = 0; nc <= nBlocks; nc++) {
			if (nc == nBlocks) {
				nread = nRem;
			}
			fread(vData, sizeof(char), n4 * nread, fi);

			for (int i = 0; i < nread; i++) {
				vLoc = vData + i * n4;
				for (int j = 0; j < n4; j++) {
					bitI = vLoc[j];
					for (int k = 0; k < 4; k++) {
						bitAB = bitI & 3;
						indTT[4 * j + k] = (bitAB >> 1) + ((bitAB & 1) << 1);
						bitI = bitI >> 2;
					}
				}
				wTT[0] = 0;
				for (int j = 0; j < ib; j++) {
					wTT[j + 1] = wTT[j] + indTT[j];
				}
				for (int j = 0; j < ib; j++) {
					wTT[ib + j + 1] = wTT[ib + j] + indTT[j];
				}
				for (int j = 0; j < ib; j++) {
					indTT[j + ib] = indTT[j];
				}
				for (int ik = 0; ik < countLen4Block[ib]; ik++) {
					ir = indexLen4Block[ib * readLM1 + ik];
					ms = ir - ib + 1;
					mB = ms / ib;
					mR = ms - mB * ib;
					sumTT = 0;

					for (int ij = 0; ij < ib; ij++) {
						sumQ = wTT[ij + mR] - wTT[ij];
						if (sumQ > sumTT)sumTT = sumQ;
					}
					sumTT += mB * wTT[ib];
					if (weightReadMax[ir] < sumTT) weightReadMax[ir] = sumTT;
				}
			}
		}

		fclose(fi); fi = nullptr;
	}

	for (int i = readLenMin; i <= readLenMax; i++) {
		printf("(2) read (%i), weight (%i)\n", i, weightReadMax[i]);
	}

	return 0;
}

int fprintString(FILE* fo, int* indTT, int n, int* n1, int* n2) {
	char sString[1000];
	int ii;
	*n1 = 0;
	*n2 = 0;
	ii = 0;
	for (int i = 0; i < n; i++) {
		if (indTT[i] == 0) {
			sString[ii] = '_'; ii++;
			continue;
		}
		if (indTT[i] == 2) {
			sString[ii] = '#'; ii++;
			(*n1)++;
		}
		else {
			sString[ii] = '@'; ii++;
			(*n2)++;
		}
	}
	sString[ii] = '\0';
	fprintf(fo, "%s\n", sString);
	return 0;
}

int fprintSeed(FILE* fo, int* indTT, int n, int mB, int mR) {
	char sString[1000];
	int ii;
	ii = 0;
	for (int m = 0; m < mB; m++) {
		for (int i = 0; i < n; i++) {
			if (indTT[i] == 0) {
				sString[ii] = '_'; ii++;
				continue;
			}
			if (indTT[i] == 2) {
				sString[ii] = '#'; ii++;
			}
			else {
				sString[ii] = '@'; ii++;
			}
		}
	}
	for (int i = 0; i < mR; i++) {
		if (indTT[i] == 0) {
			sString[ii] = '_'; ii++;
			continue;
		}
		if (indTT[i] == 2) {
			sString[ii] = '#'; ii++;
		}
		else {
			sString[ii] = '@'; ii++;
		}
	}
	sString[ii] = '\0';
	fprintf(fo, "%s\n", sString);
	return 0;
}



int findMW::reportMaxWeight() {
	int n4, rowsPerBlock, bitAB, bitI, ir;
	int ms, mB, mR, n1, n2;
	long long rowsTotal;
	int nBlocks, nRem, nread, sumTT, sumQ;
	char* vLoc;



	pfOut = (FILE**)malloc(sizeof(FILE*) * readLM1);

	for (int i = readLenMin; i <= readLenMax; i++) {
		sprintf(outputFile, "%s/tseed_mT_%i_mV_%i_readLen_%i.txt", outputFolder, mT, mV, i);
		pfOut[i] = fopen(outputFile, "w");
		fprintf(pfOut[i], "Read length: %i\n", i);
		fprintf(pfOut[i], "Mismatches (transition): %i\n", mT);
		fprintf(pfOut[i], "Mismatches (transversion): %i\n", mV);
		fprintf(pfOut[i], "Seed weight: %i", weightReadMax[i] / 2);
		if (weightReadMax[i] % 2 == 1) {
			fprintf(pfOut[i], ".5\n\n");
		}
		else {
			fprintf(pfOut[i], "\n\n");
		}
	}


	for (int ib = blockSizeMin; ib <= blockSizeMax; ib++) {

		if (countLen4Block[ib] == 0) continue;

		sprintf(inputFile, "%s/tblock_len_%i_mT_%i_mV_%i.bin", inputFolder, ib, mT, mV);
		printf("Input file (report): \"%s\"\n", inputFile);
		fi = fopen(inputFile, "rb");

		if (fi == nullptr) {
			printf("Warning: cannot open the file\n");
			continue;
		}

		n4 = ib / 4;
		if (ib % 4 > 0) n4++;

		rowsPerBlock = MEM_BLOCK_SIZE / n4;

#ifdef WIN32
		_fseeki64(fi, 0, SEEK_END);
		fileSize = _ftelli64(fi);
#else
		fseeko64(fi, 0, SEEK_END);
		fileSize = ftello64(fi);
#endif

		rowsTotal = fileSize / ((long long)n4);
		printf("Rows: %lli\n", rowsTotal);

		nBlocks = (int)(rowsTotal / ((long long)rowsPerBlock));
		nRem = int(rowsTotal - ((long long)nBlocks) * ((long long)rowsPerBlock));

		if (nRem == 0) {
			nBlocks--;
			nRem = rowsPerBlock;
		}

#ifdef WIN32
		_fseeki64(fi, 0, SEEK_SET);
#else
		fseeko64(fi, 0, SEEK_SET);
#endif

		nread = rowsPerBlock;

		for (int nc = 0; nc <= nBlocks; nc++) {
			if (nc == nBlocks) {
				nread = nRem;
			}
			fread(vData, sizeof(char), n4 * nread, fi);

			for (int i = 0; i < nread; i++) {
				vLoc = vData + i * n4;
				for (int j = 0; j < n4; j++) {
					bitI = vLoc[j];
					for (int k = 0; k < 4; k++) {
						bitAB = bitI & 3;
						indTT[4 * j + k] = (bitAB >> 1) + ((bitAB & 1) << 1);
						bitI = bitI >> 2;
					}
				}
				wTT[0] = 0;
				for (int j = 0; j < ib; j++) {
					wTT[j + 1] = wTT[j] + indTT[j];
				}
				for (int j = 0; j < ib; j++) {
					wTT[ib + j + 1] = wTT[ib + j] + indTT[j];
				}
				for (int j = 0; j < ib; j++) {
					indTT[j + ib] = indTT[j];
				}
				for (int ik = 0; ik < countLen4Block[ib]; ik++) {
					ir = indexLen4Block[ib * readLM1 + ik];
					ms = ir - ib + 1;
					mB = ms / ib;
					mR = ms - mB * ib;
					sumTT = weightReadMax[ir] - mB * wTT[ib];
					if (sumTT > wTT[ib])continue;

					for (int ij = 0; ij < ib; ij++) {
						sumQ = wTT[ij + mR] - wTT[ij];
						//if (indTT[ij] == 0)continue;
						if (sumQ < sumTT)continue;

						fprintf(pfOut[ir], "Block size: %i\n", ib);
						fprintf(pfOut[ir], "Number of whole blocks: %i\n", mB);
						fprintf(pfOut[ir], "Remainder: %i\n", mR);
						fprintString(pfOut[ir], indTT + ij, ib, &n1, &n2);
						fprintf(pfOut[ir], "Number of # in the block: %i\n", n1);
						fprintf(pfOut[ir], "Number of @ in the block: %i\n", n2);
						fprintSeed(pfOut[ir], indTT + ij, ib, mB, mR);
						fprintf(pfOut[ir], "\n");
					}
				}
			}
		}

		fclose(fi); fi = nullptr;
	}

	return 0;
}

int main(int argc, char** argv) {
	findMW* fmw;
	int ires;

	fmw = new findMW();

	ires = fmw->startProcessing(argc, argv);

	delete fmw; fmw = nullptr;

	return ires;
}
