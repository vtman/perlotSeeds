#include <fstream>

const int MEM_BLOCK_SIZE = 10000000;

int printInfoParameters() {
	printf("Input parameters:\n");
	printf("1) Input (binary) folder\n");
	printf("2) Output (text) folder\n");
	printf("3) Block size\n");
	printf("4) Number of transitions\n");
	printf("5) Number of transversions\n\n");
	return 0;
}

int main(int argc, char** argv) {
	int mT, mV, blockSize, nbr, rowsPerBlock, nB, nR, iValue;
	char inputFolder[1000], outputFolder[1000], inputFile[1000], outputFile[1000], sLine[1000];
	unsigned char *outVector, *uRes;
	long long fileSize, numRows;
	FILE* fi, * fo;

	if (argc != 6) {
		printInfoParameters();
		return -1;
	}

	sprintf(inputFolder, "%s", argv[1]);
	sprintf(outputFolder, "%s", argv[2]);
	blockSize = atoi(argv[3]);
	mT = atoi(argv[4]);
	mV = atoi(argv[5]);

	printf("\nInput (binary) folder: %s\n", inputFolder);
	printf("Output (text) folder: %s\n", outputFolder);
	printf("Block size: %i\n", blockSize);
	printf("Number of transitions: %i\n", mT);
	printf("Number of transversions: %i\n\n", mV);

	sprintf(inputFile, "%s/tblock_len_%i_mT_%i_mV_%i.bin", inputFolder, blockSize, mT, mV);
	printf("Input file: \"%s\"\n", inputFile);

	nbr = blockSize / 4;
	if ((blockSize & 3) > 0) nbr++;

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open the input file\n");
		return -2;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);

	numRows = fileSize / ((long long)nbr);
	rowsPerBlock = MEM_BLOCK_SIZE / nbr;

	printf("Number of rows: %lli\n", numRows);

	nB = numRows / ((long long)rowsPerBlock);
	nR = (int)(numRows - (long long)nB * (long long)rowsPerBlock);

	sprintf(outputFile, "%s/tblock_len_%i_mT_%i_mV_%i.txt", outputFolder, blockSize, mT, mV);
	printf("Output file: \"%s\"\n", outputFile);

	fo = fopen(outputFile, "w");
	if (fo == nullptr) {
		printf("Error: cannot open the output file\n");
		fclose(fi); fi = nullptr;
		return -3;
	}


	outVector = (unsigned char*)malloc(sizeof(char) * MEM_BLOCK_SIZE);

	sLine[blockSize] = '\0';

	for (int ib = 0; ib < nB; ib++) {
		fread(outVector, sizeof(char), nbr * rowsPerBlock, fi);
		for (int i = 0; i < rowsPerBlock; i++) {
			uRes = outVector + i * nbr;
			for (int j = 0; j < blockSize; j++) {
				iValue = 3 & (uRes[j >> 2] >> (2 * (j & 3)));
				if (iValue == 1) {
					sLine[j] = '#';
				}
				else if (iValue == 2) {
					sLine[j] = '@';
				}
				else {
					sLine[j] = '_';
				}
			}
			fprintf(fo, "%s\n", sLine);
		}
	}

	if (nR > 0) {
		fread(outVector, sizeof(char), nbr * nR, fi);
		for (int i = 0; i < nR; i++) {
			uRes = outVector + i * nbr;
			for (int j = 0; j < blockSize; j++) {
				iValue = 3 & (uRes[j >> 2] >> (2 * (j & 3)));
				if (iValue == 1) {
					sLine[j] = '#';
				}
				else if (iValue == 2) {
					sLine[j] = '@';
				}
				else {
					sLine[j] = '_';
				}
			}
			fprintf(fo, "%s\n", sLine);
		}
	}

	fclose(fo); fo = nullptr;
	fclose(fi); fi = nullptr;

	free(outVector); outVector = nullptr;

	return 0;
}