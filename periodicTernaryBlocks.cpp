//E:\Temp2\perlotSeeds\binaryBlocks E:\Temp2\perlotSeeds\ternaryBlocks 30 2 3 0

#include <fstream>

#include <chrono>
#include <iostream>

#include <smmintrin.h>
#include <emmintrin.h>

#include <string.h>

#include <omp.h>

const int MAX_ELEMENTS = 1000000;

class blockData {
public:
	blockData(int nT_in, int nCheck_in);
	~blockData();

	int processLine(int lineIndex);

	int nT, nCheck, nResult;
	bool isSmallSize;
	__m128i* v1, * v2, * vResult;

};

class TBSeed {
public:
	TBSeed();
	~TBSeed();
	int startProcessing(int narg, char** varg);
	int printUsage();
	int printInputInfo();
	int getWeight(int iw, int* wRes);

	int printM(__m128i* vec);
	int printTT(__m128i mA, __m128i mB);
	int printTTbin(__m128i mA, __m128i mB);
	int getDataTransition();
	int getDataTot();

	int freeDataTransition();
	int freeDataTotal();

	blockData** pbd;

	__m128i** pDataTransition, ** pDataTotal;
	int nbTransition, nbc, nbTotal, nRemTransition, nRemTotal;

	FILE* fi, * fo;
	char inputFolder[1000], inputFile[1000], outputFolder[1000], outputFile[1000];
	char sString[1000], oString[1000];
	unsigned char* vD1, * vD2, * vD3;
	int iOut[1000];
	unsigned char uOut[1000];
	int mT, mV, mTot, blockSize, iLevel, iLevelT, iLevelTot, L8;
	int weightMin, weightMax, nThreads, nT, nCheck, nResult;
	bool resST, isFileClosed;

	__m128i* vmT, * vmCheck, * vmResult, mOne[128], mZero;

	long long nRowTransition, nRowTotal;
};

blockData::blockData(int nT_in, int nCheck_in) {
	nT = nT_in;
	nCheck = nCheck_in;

	v1 = (__m128i*)malloc(sizeof(__m128i) * nT);
	v2 = (__m128i*)malloc(sizeof(__m128i) * nCheck);
	vResult = (__m128i*)malloc(2 * sizeof(__m128i) * MAX_ELEMENTS);

	nResult = 0;
}

blockData::~blockData() {
	if (v1 != nullptr) { free(v1); v1 = nullptr; }
	if (v2 != nullptr) { free(v2); v2 = nullptr; }
	if (vResult != nullptr) { free(vResult); vResult = nullptr; }
}

int blockData::processLine(int lineIndex) {
	__m128i m1, m2, m3;
	m1 = v1[lineIndex];

	if (!isSmallSize) {
		for (int i = 0; i < nCheck; i++) {
			m2 = _mm_or_si128(m1, v2[i]);
			m3 = _mm_xor_si128(m1, m2);
			if (_mm_extract_epi64(m3, 0) == 0) {
				vResult[2 * nResult] = m1;
				vResult[2 * nResult + 1] = v2[i];
				nResult++;
				if (nResult > MAX_ELEMENTS) {
					printf("Error: small size storage\n");
					isSmallSize = true;
					break;
				}
			}
		}
	}
	return 0;
}


TBSeed::TBSeed() {
	int a[4], i1, i2;

	pDataTransition = nullptr;
	pDataTotal = nullptr;
	nbTotal = 0;
	nbTransition = 0;

	vmCheck = nullptr;
	vmResult = nullptr;
	vmT = nullptr;
	vD1 = nullptr;
	vD2 = nullptr;
	vD3 = nullptr;
	fi = nullptr;
	fo = nullptr;
	resST = false;

	isFileClosed = true;

	pbd = nullptr;
	nThreads = 0;
	
	for (int i = 0; i < 128; i++) {
		i1 = i / 32;
		i2 = i % 32;
		a[0] = 0;
		a[1] = 0;
		a[2] = 0;
		a[3] = 0;

		a[i1] = 1 << i2;
		mOne[i] = _mm_set_epi32(a[3], a[2], a[1], a[0]);
	}
	mZero = _mm_set1_epi32(0);
}

int TBSeed::freeDataTransition() {
	if (pDataTransition != nullptr) {
		for (int i = 0; i < nbTransition; i++) {
			if (pDataTransition[i] != nullptr) {
				free(pDataTransition[i]); pDataTransition[i] = nullptr;
			}
		}
		free(pDataTransition); pDataTransition = nullptr;
	}
	return 0;
}

int TBSeed::freeDataTotal() {
	if (pDataTotal != nullptr) {
		for (int i = 0; i < nbTotal; i++) {
			if (pDataTotal[i] != nullptr) {
				free(pDataTotal[i]); pDataTotal[i] = nullptr;
			}
		}
		free(pDataTotal); pDataTotal = nullptr;
	}
	return 0;
}

TBSeed::~TBSeed() {
	freeDataTransition();
	freeDataTotal();

	if (vmT != nullptr) { free(vmT); vmT = nullptr; }
	if (vmResult != nullptr) { free(vmResult); vmResult = nullptr; }
	if (vmCheck != nullptr) { free(vmCheck); vmCheck = nullptr; }
	if (pbd != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			if (pbd[i] != nullptr) {
				delete pbd[i]; pbd[i] = nullptr;
			}
		}
		delete pbd; pbd = nullptr;
	}
	if (vD1 != nullptr) { free(vD1); vD1 = nullptr; }
	if (vD2 != nullptr) { free(vD2); vD2 = nullptr; }
	if (vD3 != nullptr) { free(vD3); vD3 = nullptr; }

	if (fi != nullptr) { fclose(fi); fi = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr; }
}

int TBSeed::getWeight(int iw, int* wRes) {
	int nt, isq;
	if (iw == 1) {
		*wRes = blockSize - 1;
		return 0;
	}
	sprintf(inputFile, "%s/block_m_%i_len_%i_level_0.txt", inputFolder, iw, blockSize);
	fi = fopen(inputFile, "r");
	if (fi == nullptr) {
		printf("Error: cannot read from file %s\n", inputFile);
		return -1;
	}
	fscanf(fi, "%s", sString);
	fclose(fi); fi = nullptr;

	nt = 0;
	isq = strlen(sString);
	printf("%s|\n", sString);
	for (int i = 0; i < isq; i++) {
		if (sString[i] == '1') nt++;
	}
	printf("Weight: %i\n", nt);
	*wRes = nt;
	return 0;
}


int TBSeed::getDataTransition() {
	long long fSize;
	int nread, iA, aa[4];
	char* vData;
	__m128i mA;

	if (iLevelT == 0 && mT == 1) {
		nbTransition = 1;
		nRemTransition = 1;
		nRowTransition = 1;

		pDataTransition = (__m128i**)malloc(sizeof(__m128i*) * nbTransition);

		pDataTransition[0] = (__m128i*)malloc(sizeof(__m128i) * 1);
		aa[0] = 0;
		aa[1] = 0;
		aa[2] = 0;
		aa[3] = 0;

		for (int i = 1; i < blockSize; i++) {
			aa[i / 32] |= 1 << (i & 31);
		}

		pDataTransition[0][0] = _mm_set_epi32(aa[3], aa[2], aa[1], aa[0]);

		printf("nRow (transition): %lli\n\n", nRowTransition);

		return 0;
	}

	sprintf(inputFile, "%s/block_m_%i_len_%i_level_%i.bin", inputFolder, mT, blockSize, iLevelT);
	printf("Transition data: \"%s\"\n", inputFile);
	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot read from file %s\n", inputFile);
		return -1;
	}


#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fSize = _ftelli64(fi);
#else
	fseeko64(fi, 0, SEEK_END);
	fSize = ftello64(fi);
#endif

	nRowTransition = fSize / (long long)(L8);
	rewind(fi);

	nbTransition = nRowTransition / nbc;
	nRemTransition = nRowTransition % nbc;
	if (nRemTransition > 0) {
		nbTransition++;
	}
	else {
		nRemTransition = nbc;
	}

	pDataTransition = (__m128i**)malloc(sizeof(__m128i*) * nbTransition);

	nread = nbc;

	vData = (char*)malloc(sizeof(char) * nbc * L8);

	for (int i = 0; i < nbTransition; i++) {
		if (i < nbTransition - 1) {
			pDataTransition[i] = (__m128i*)malloc(sizeof(__m128i) * nbc);
		}
		else {
			pDataTransition[i] = (__m128i*)malloc(sizeof(__m128i) * nRemTransition);
			nread = nRemTransition;
		}
		fread(vData, sizeof(char), L8 * nread, fi);

		for (int j = 0; j < nread; j++) {
			mA = mZero;
			for (int k = 0; k < L8; k++) {
				iA = vData[j * L8 + k];
				for (int ii = 0; ii < 8; ii++) {
					if ((iA & 1) == 1) {
						mA = _mm_or_si128(mA, mOne[k * 8 + ii]);
					}
					iA = iA >> 1;
				}
			}
			pDataTransition[i][j] = mA;
		}
	}

	free(vData); vData = nullptr;

	fclose(fi); fi = nullptr;

	printf("nRow (transition): %lli\n\n", nRowTransition);
	return 0;
}

bool areSame(__m128i m1, __m128i m2) {
	__m128i m3, m4, m5;
	m3 = _mm_xor_si128(m1, m2);
	m4 = _mm_bsrli_si128(m3, 8);
	m5 = _mm_or_si128(m3, m4);
	if (_mm_extract_epi64(m5, 0) == 0) return true;
	return false;
}


int TBSeed::getDataTot() {
	long long fSize, nQ;
	int nread, nC, nChunk, nRem, nCount, kk, iValue, nCount2;
	bool t;
	char* vData;
	__m128i mA, mM[200];

	sprintf(inputFile, "%s/block_m_%i_len_%i_level_%i.bin", inputFolder, mTot, blockSize, iLevelTot);
	printf("Total data: \"%s\"\n", inputFile);
	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot read from file %s\n", inputFile);
		return -1;
	}

#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fSize = _ftelli64(fi);
#else
	fseeko64(fi, 0, SEEK_END);
	fSize = ftello64(fi);
#endif

	nRowTotal = fSize / (long long)(L8);

	printf("Total rows: %lli\n", nRowTotal);
	rewind(fi);

	nChunk = nRowTotal / nbc;
	nRem = nRowTotal % nbc;
	if (nRem > 0) {
		nChunk++;
	}
	else {
		nRem = nbc;
	}

	nQ = (long long)(2 * blockSize) * nRowTotal;
	nC = nQ / nbc;
	if (nQ % nbc > 0) nC++;

	pDataTotal = (__m128i**)malloc(sizeof(__m128i*) * nC);
	for (int i = 0; i < nC; i++) {
		pDataTotal[i] = nullptr;
	}

	nread = nbc;

	vData = (char*)malloc(sizeof(char) * nbc * L8);

	nbTotal = 0;
	nRemTotal = 0;

	pDataTotal[0] = (__m128i*)malloc(sizeof(__m128i) * nbc);

	for (int i = 0; i < nChunk; i++) {
		if (i == nChunk - 1) {
			nread = nRem;
		}
		fread(vData, sizeof(char), L8 * nread, fi);

		for (int j = 0; j < nread; j++) {

			nCount = 0;
			for (int k = 0; k < blockSize; k++) {
				mA = mZero;
				iValue = vData[j * L8 + (k >> 3)];
				if (((iValue >> (k & 7)) & 1) == 1)continue;
				for (int m = 0; m < blockSize; m++) {
					kk = k + m;
					if (kk >= blockSize) kk -= blockSize;
					iValue = vData[j * L8 + (kk >> 3)];
					if ((1 & (iValue >> (kk & 7))) == 1) {
						mA = _mm_or_si128(mA, mOne[m]);
					}
				}
				mM[nCount] = mA;
				nCount++;
			}

			for (int k = 0; k < blockSize; k++) {
				mA = mZero;
				iValue = vData[j * L8 + (k >> 3)];
				if (((iValue >> (k & 7)) & 1) == 1)continue;
				for (int m = 0; m < blockSize; m++) {
					kk = k - m;
					if (kk < 0) kk += blockSize;
					iValue = vData[j * L8 + (kk >> 3)];
					if ((1 & (iValue >> (kk & 7))) == 1) {
						mA = _mm_or_si128(mA, mOne[m]);
					}
				}
				mM[nCount] = mA;
				nCount++;
			}

			nCount2 = 1;
			for (int k = 1; k < nCount; k++) {
				t = false;
				for (int ij = 0; ij < nCount2; ij++) {
					if (areSame(mM[ij], mM[k])) {
						t = true;
						break;
					}
				}
				if (t) {
					continue;
				}
				mM[nCount2] = mM[k];
				nCount2++;
			}

			for (int k = 0; k < nCount2; k++) {
				pDataTotal[nbTotal][nRemTotal] = mM[k];
				nRemTotal++;
				if (nRemTotal == nbc) {
					nRemTotal = 0;
					nbTotal++;
					pDataTotal[nbTotal] = (__m128i*)malloc(sizeof(__m128i) * nbc);
				}
			}

		}
	}

	free(vData); vData = nullptr;

	fclose(fi); fi = nullptr;

	if (nRemTotal == 0) {
		nRemTotal = nbc;
	}
	else {
		nbTotal++;
	}

	nRowTotal = (long long)nbc * (long long)(nbTotal - 1) + (long long)nRemTotal;

	printf("nRow (total): %lli\n\n", nRowTotal);
	return 0;
}

int TBSeed::startProcessing(int narg, char** varg) {
	bool isFileOpen;
	int nr1, nr2;
	__m128i* pv1, * pv2, mA, mB, mC;

	if (narg != 7) {
		printf("Error: wrong number of arguments\n");
		printUsage();
		return -1;
	}
	sprintf(inputFolder, "%s", varg[1]);
	sprintf(outputFolder, "%s", varg[2]);
	blockSize = atoi(varg[3]);
	mT = atoi(varg[4]);
	mV = atoi(varg[5]);
	iLevel = atoi(varg[6]);
	mTot = mT + mV;

	if (printInputInfo() != 0) return -2;

	L8 = blockSize / 8;
	if (blockSize % 8 > 0)L8++;

	isFileOpen = false;

	nbc = MAX_ELEMENTS;


	for (int iL = 0; iL <= iLevel; iL++) {
		iLevelT = iL;
		iLevelTot = iLevel - iL;
		printf("Level (transition): %i, level (total): %i\n\n", iLevelT, iLevelTot);

		if (getDataTransition() != 0) return -3;
		if (getDataTot() != 0) return -4;

		for (int ib1 = 0; ib1 < nbTransition; ib1++) {
			printf("Block #%i of %i\n", ib1, nbTransition);
			if (ib1 < nbTransition - 1) {
				nr1 = nbc;
			}
			else {
				nr1 = nRemTransition;
			}
			pv1 = pDataTransition[ib1];

			for (int i = 0; i < nr1; i++) {
				mA = pv1[i];

				for (int ib2 = 0; ib2 < nbTotal; ib2++) {
					if (ib2 < nbTotal - 1) {
						nr2 = nbc;
					}
					else {
						nr2 = nRemTotal;
					}
					pv2 = pDataTotal[ib2];

					for (int j = 0; j < nr2; j++) {
						mB = pv2[j];
						mC = _mm_or_si128(mA, mB);
						if (!areSame(mA, mC))continue;

						if (!isFileOpen) {
							sprintf(outputFile, "%s/tblock_len_%i_mT_%i_mV_%i.bin", outputFolder, blockSize, mT, mV);
							fo = fopen(outputFile, "wb");
							isFileOpen = true;
						}
						printTTbin(mA, mB);
					}
				}
			}
		}

		freeDataTransition();
		freeDataTotal();
	}

	if (!isFileOpen) {
		sprintf(outputFile, "%s/tblock_len_%i_mT_%i_mV_%i_no.txt", outputFolder, blockSize, mT, mV);
		fo = fopen(outputFile, "wb");
	}

	fclose(fo); fo = nullptr;
	return 0;
}

int TBSeed::printInputInfo() {
	printf("\nInput folder: %s\n", inputFolder);
	printf("Output folder: %s\n", outputFolder);
	printf("Block size: %i\n", blockSize);
	printf("Number of mismatches (transition): %i\n", mT);
	printf("Number of mismatches (transversion): %i\n", mV);
	printf("Level: %i\n\n", iLevel);
	return 0;
}

int TBSeed::printUsage() {
	printf("Arguments:\n");
	printf("\t1) Input folder\n");
	printf("\t2) Output folder\n");
	printf("\t3) Block size\n");
	printf("\t4) Number of mismatches (transition)\n");
	printf("\t5) Number of mismatches (transversion)\n");
	printf("\t6) Level\n\n");

	return 0;
}

int TBSeed::printM(__m128i* vec) {
	int a[4], b[4], i1, i2;
	int bitA, bitB;

	_mm_store_si128((__m128i*)a, vec[0]);
	_mm_store_si128((__m128i*)b, vec[1]);

	for (int i = 0; i < blockSize; i++) {
		i1 = i / 32;
		i2 = i % 32;
		bitA = ((a[i1] >> i2) & 1);
		bitB = ((b[i1] >> i2) & 1);
		if (bitB == 1) {
			fprintf(fo, "@");
		}
		else {
			if (bitA == 1) {
				fprintf(fo, "#");
			}
			else {
				fprintf(fo, "_");
			}
		}
	}
	fprintf(fo, "\n");

	return 0;
}

int TBSeed::printTT(__m128i mA, __m128i mB) {
	int a[4], b[4], i1, i2;
	int bitA, bitB;

	_mm_store_si128((__m128i*)a, mA);
	_mm_store_si128((__m128i*)b, mB);

	for (int i = 0; i < blockSize; i++) {
		i1 = i / 32;
		i2 = i % 32;
		bitA = ((a[i1] >> i2) & 1);
		bitB = ((b[i1] >> i2) & 1);
		if (bitB == 1) {
			fprintf(fo, "#");
		}
		else {
			if (bitA == 1) {
				fprintf(fo, "@");
			}
			else {
				fprintf(fo, "_");
			}
		}
	}
	fprintf(fo, "\n");

	return 0;
}

int TBSeed::printTTbin(__m128i mA, __m128i mB) {
	int a[4], b[4], i1, i2, i3, i4;
	int bitA, bitB, b4, bitC;
	unsigned char u[100];
	int aa[100];

	b4 = blockSize >> 2;
	if ((blockSize & 3) > 0)b4++;

	_mm_store_si128((__m128i*)a, mA);
	_mm_store_si128((__m128i*)b, mB);

	for (int i = 0; i < b4; i++) {
		aa[i] = 0;
	}

	for (int i = 0; i < blockSize; i++) {
		i1 = i >> 5;
		i2 = i & 31;
		i3 = i >> 2;
		i4 = (i & 3) << 1;
		bitA = ((a[i1] >> i2) & 1);
		bitB = ((b[i1] >> i2) & 1);
		bitC = ((bitA ^ bitB) << 1) | bitB;

		aa[i3] |= (bitC << i4);
	}

	for (int i = 0; i < b4; i++) {
		u[i] = (unsigned char)aa[i];
	}

	fwrite(u, sizeof(char), b4, fo);

	return 0;
}


int main(int argc, char** argv) {
	int ires;
	TBSeed* ts;

	auto start = std::chrono::steady_clock::now();

	ts = new TBSeed();
	ires = ts->startProcessing(argc, argv);
	delete ts; ts = nullptr;

	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("It took me %lld ms.\n\n", elapsed.count());

	return ires;
}
