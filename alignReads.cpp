#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <chrono>

#include <omp.h>

#include <emmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>

//#define WIN32 1

//#define CASE_T0V1 1
//#define CASE_T0V2 1
//#define CASE_T0V3 1
//#define CASE_T1V0 1
//#define CASE_T1V1 1
//#define CASE_T1V2 1
//#define CASE_T1V3 1
//#define CASE_T2V0 1
//#define CASE_T2V1 1
//#define CASE_T2V2 1
//#define CASE_T3V0 1
//#define CASE_T3V1 1
//#define CASE_T3V2 1
//#define CASE_T4V0 1
//#define CASE_T4V1 1
//#define CASE_T5V0 1
//#define CASE_T5V1 1
//#define CASE_T6V0 1
//#define CASE_C16 1
//#define CASE_C20 1
//#define CASE_C24 1
//#define CASE_C28 1
//#define CASE_C32 1
//#define CASE_C36 1
#define CASE_C40 1
//#define CASE_C44 1
//#define CASE_C48 1
//#define CASE_C52 1
//#define CASE_C56 1
//#define CASE_C60 1
//#define CASE_C64 1

const int MAX_READ_COUNT = 10000;
const int PRINT_READ_SHIFT = 40;
const int PRINT_EXTRA_REF = 10;

const int MAX_RECORDS = 500000;

class patternReadSeed {
public:
	patternReadSeed(char* vSeed_in, int nLetters_in, int nIndexLevel_in);
	~patternReadSeed();
	int copyData(__m128i* v_in, int iShift_in);
	int formMaskRef();
	int formSignature();
	int findWeight();
	int findFix();
	int findLibraryList(char* vc, int indFirst, int indLast, int* indListFirst, int* indListLast);
	int comparePatternLibrary(char* vcPattern, char* vcRecord);

	int nSearchBits, nSearchBytes, nSearchRem, nSearchWhole;

	__m128i* vm, * vMask, * vRef, * vPattern, * vPatternShift, * vOut;
	int nLen, nLen128, iShift, nDoubleWeight, recordSize;
	int vOne[32];
	int nLetters, nIndexLevel, nPattern;
	int indFile, indIndex;
	int uFileMask, shiftBitFile;
	bool isOk;
	char* vQ, vSeed[1000];
	__m128i m1, m2, m3, m4, m5, mZero, mTail;
	__m128i mAC, mAG, mOut;
};

class blockData {
public:
	blockData();
	~blockData();

	int processReads(int indRelativeStart_in, int nReadsToProcess_in);
	int formPatternStartPositions();
	int mergeCandiateListSpaced(int ind);
	int measureSimilarity(unsigned int upos, __m128i* vR, int* countT, int* countV);
	int findPairedList();
	int formStringPosition(unsigned int upos, char* vC);

	int* voCountCandidatesLoc, * voCountSolutionsLoc, * voBestLoc;
	int* indPatternStartPosition, * vIndexCurrent;
	int* vCountCandidate, vCountMergedCandidate[4];
	int** pvLibraryStart;

	int minPairDistance, maxPairDistance, nPatterns, recordSize;
	int nChromosomes, nCountPair03, nCountPair12;
	int indRelativeStart, nReadsToProcess, blockMsize, readMsize;

	unsigned int** pvCandidatePosition, ** pvMergedCandidatePosition;
	unsigned int* vPairPosition03, * vPairPosition12;

	char** pvLibraryData;

	__m128i* vmReadLeft, * vmReadRight, * vmReadLoc, * mSim, * vmRef;

	long long indReadStart;
	long long* uChromStart;

	patternReadSeed** pPRS;

	FILE* fo;
};

class contP {
public:
	contP();
	~contP();

	int loadLibrary();

	int printInfo();
	int startProcessing(int nArg, char **vArg);
	int setParameters();
	int getSeedInfo();
	int getReadInfo();

	int printReference(unsigned int uPos);
	int getReferenceData();
	int printLibraryRecord(char* vc, int ind);
	int formPatternStartPositions();

	int getInfoReferenceData();

	blockData** pBD;

	char** pvLibraryData;
	char libraryTemplate[1000], printReadTemplateStart[1000], printReadTemplateGap[1000];
	char printRefTemplateStart[1000], printRefTemplateGap[1000];
	char inputPairLeft[1000], inputPairRight[1000], inputInfoReference[1000], inputReference[1000], inputLibrary[1000];
	char inputFile[1000], outputFolder[1000], outputFile[1000], sSeed[1000], sTemp[1000];
	char stringRef[1000], stringRead[1000];

	unsigned int** pvPositions;

	int** pvLibraryStart;
	int* voCountCandidates, * voCountSolutions, * voBest;
	int* vCountPositions, * vnLibraryRecords;
	int maxPairDistance, minPairDistance, nPatterns, nThreads, nChromosomes, recordSize;
	int indRelativeStart, nReadsToProcess, countRR;
	int seedLength, seedDoubleWeight, seedNLetters, seedIndexLevel;
	int readLengthLeft, readLengthRight, readBinarySize, readBlockSize;
	int readMsize, blockMsize, readFirst, readStep;
	int nLibraryEntries, nLibraryFiles, nLenRef;

	FILE* fiLeft, * fiRight, * fiACGT, * fo, * fi, * foStat;

	long long countReadLeft, countReadRight, fileSizeLeft, fileSizeRight, indReadStart;
	long long* uChromStart;

	__m128i* vmRef;
};


blockData::blockData() {
	nPatterns = 0;
	fo = nullptr;
	vmReadLeft = nullptr;
	vmReadRight = nullptr;

	vPairPosition03 = nullptr;
	vPairPosition12 = nullptr;
	vIndexCurrent = nullptr;
	vCountCandidate = nullptr;
	pvCandidatePosition = nullptr;
	pvMergedCandidatePosition = nullptr;
	indPatternStartPosition = nullptr;

	uChromStart = nullptr;

	pPRS = nullptr;

	mSim = nullptr;

	voBestLoc = nullptr;
	voCountCandidatesLoc = nullptr;
	voCountSolutionsLoc = nullptr;
}

blockData::~blockData() {
	if (uChromStart != nullptr) { free(uChromStart); uChromStart = nullptr; }
	if(voBestLoc != nullptr){free(voBestLoc); voBestLoc = nullptr;}
	if(voCountCandidatesLoc != nullptr){free(voCountCandidatesLoc); voCountCandidatesLoc = nullptr;}
	if(voCountSolutionsLoc != nullptr){free(voCountSolutionsLoc); voCountSolutionsLoc = nullptr;}

	if (fo != nullptr) { fclose(fo); fo = nullptr; }
#ifdef WIN32
	if (vmReadLeft != nullptr) { _aligned_free(vmReadLeft); vmReadLeft = nullptr; }
	if (vmReadRight != nullptr) { _aligned_free(vmReadRight); vmReadRight = nullptr; }
	if (mSim != nullptr) { _aligned_free(mSim); mSim = nullptr; }
#else
	if (vmReadLeft != nullptr) { free(vmReadLeft); vmReadLeft = nullptr; }
	if (vmReadRight != nullptr) { free(vmReadRight); vmReadRight = nullptr; }
	if (mSim != nullptr) { free(mSim); mSim = nullptr; }
#endif

	if (vIndexCurrent != nullptr) { free(vIndexCurrent); vIndexCurrent = nullptr; }
	if (indPatternStartPosition != nullptr) {
		free(indPatternStartPosition); indPatternStartPosition = nullptr;
	}

	if (pvCandidatePosition != nullptr) {
		for (int i = 0; i < 4 * nPatterns; i++) {
			if (pvCandidatePosition[i] != nullptr) {
				free(pvCandidatePosition[i]); pvCandidatePosition[i] = nullptr;
			}
		}
		free(pvCandidatePosition); pvCandidatePosition = nullptr;
	}
	if (vPairPosition03 != nullptr) {
		free(vPairPosition03); vPairPosition03 = nullptr;
	}
	if (vPairPosition12 != nullptr) {
		free(vPairPosition12); vPairPosition12 = nullptr;
	}

	if (pvMergedCandidatePosition != nullptr) {
		for (int i = 0; i < 4; i++) {
			if (pvMergedCandidatePosition[i] != nullptr) {
				free(pvMergedCandidatePosition[i]); pvMergedCandidatePosition[i] = nullptr;
			}
		}
		free(pvMergedCandidatePosition); pvMergedCandidatePosition = nullptr;
	}
	if (vCountCandidate != nullptr) { free(vCountCandidate); vCountCandidate = nullptr; }

	if (pPRS != nullptr) {
		for (int i = 0; i < nPatterns; i++) {
			if (pPRS[i] != nullptr) {
				delete pPRS[i]; pPRS[i] = nullptr;
			}
		}
		free(pPRS); pPRS = nullptr;
	}
}


int formSignatureM(__m128i* m, __m128i* vPattern) {

#ifdef CASE_T0V1
	//T0V1
	__m128i c, t, s, mAC, mAG, mMask, res[2];

	c = _mm_set1_epi32(0xbfbfbfbf);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 25);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x20000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 15);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000010);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000008);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 27);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x1fbfbfbf);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00000004);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 4);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x00000002);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 13);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x00000001);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 22);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x1fffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 29));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_T0V2
	//T0V2
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xf4efa77d);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x40000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 29);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x20110000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x1600a000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00001e9d);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00800000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 18);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x01420000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00001fff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 13));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_T0V3
	//T0V3
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xc8f591eb);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00200040);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x32086410);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00140028);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00000002);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00010000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00000003);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 2));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

#endif

#ifdef CASE_T1V0
	//T1V0
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xbfbfbfbf);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 25);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x20000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 15);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000010);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000008);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 27);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x1fbfbfbf);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00000004);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 4);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x00000002);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 13);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x00000001);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 22);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x1fffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 29));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
	vPattern[1] = _mm_set1_epi32(0);


	c = _mm_set1_epi32(0x00000040);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x40000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 30);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00400000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 21);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00004000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 12);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x40000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 27);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00400000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00004000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000040);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 26), 12));
	vPattern[1] = _mm_or_si128(vPattern[1], _mm_srli_epi32(res[0], 6));
#endif

#ifdef CASE_T1V1
	//T1V1
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x9f59f59f);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80180000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 14);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x14014000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x60040000);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);

	c = _mm_set1_epi32(0x000019f5);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x01820000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 14);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x00000001);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 1);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00001fff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 13));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	c = _mm_set1_epi32(0x00000a40);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x40000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 23);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00a40000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0a008000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 13);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00402400);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0000000a);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 26), 8));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 6), 12));
#endif

#ifdef CASE_T1V2
	//T1V2

	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xc8f591eb);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00200040);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x32086410);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00140028);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00000002);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00010000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00000003);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 2));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);


	c = _mm_set1_epi32(0x00000000);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x14000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 26);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00002800);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 10);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x05000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 20);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000a00);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 4), 8));
#endif

#ifdef CASE_T1V3
	__m128i c, t, s, mAC, mAG, mMask, res[2];

	c = _mm_set1_epi32(0x00219143);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0xc8800000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 20);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x32002000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 11);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00204440);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00080010);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);

	mMask = _mm_set_epi32(0, 0, 0, 0x003fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 22));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));

	c = _mm_set1_epi32(0x000000a8);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00540000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00100000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 19);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00050000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0000002a);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 12), 4));
#endif

#ifdef CASE_T2V0
	//T2V0
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xf4efa77d);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x40000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 29);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x20110000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x1600a000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00001e9d);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00800000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 18);
	res[1] = _mm_or_si128(res[1], s);
	c = _mm_set1_epi32(0x01420000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00001fff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 13));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	c = _mm_set1_epi32(0x00005882);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x02000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 20);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x09100000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 18);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x08284000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 11);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00040022);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000140);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 26), 8));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 6), 12));
#endif

#ifdef CASE_T2V1
	//T2V1

	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xc8f591eb);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00200040);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x32086410);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00140028);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00000002);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00010000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00000003);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 2));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	c = _mm_set1_epi32(0x000a6c14);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x30000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 28);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x06000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 10);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x04820900);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x09001200);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00008005);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 4), 8));
#endif

#ifdef CASE_T2V2
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x0060358b);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 29);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x31000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00d20c00);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00042002);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000004);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x00ffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 24));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));

	c = _mm_set1_epi32(0x000cc074);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x0e800000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 15);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0000c000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 14);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00010088);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00001110);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 16), 4));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 16), 8));

#endif

#ifdef CASE_T3V0
	//T3V0
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0xc8f591eb);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x00200040);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x32086410);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00140028);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x00000002);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00010000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 16);
	res[1] = _mm_or_si128(res[1], s);

	mAC = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4));
	mAG = _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	mMask = _mm_set_epi32(0, 0, 0, 0x00000003);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 2));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);


	c = _mm_set1_epi32(0x070a6e14);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x30000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 28);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0d001200);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00808100);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00420885);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 4), 8));
#endif

#ifdef CASE_T3V1
	//T3V1
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x00110c43);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x04000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 24);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0xc0200000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 16);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x31044000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 11);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00080310);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00000004);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 10);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x001fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 21));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));



	c = _mm_set1_epi32(0x33ccf33c);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0xc0330c00);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 10);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0cc03000);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x000000cf);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x00000003);
	t = _mm_and_si128(m[2], c);
	s = _mm_slli_epi32(t, 4);
	res[1] = _mm_or_si128(res[1], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[1] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 10), 4));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 22), 8));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[1], 10), 8));
#endif

#ifdef CASE_T3V2
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x00001219);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x08000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 26);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x40640000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 16);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x06410000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00002100);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 2);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000090);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x0003ffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 18));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));



	c = _mm_set1_epi32(0x009ac526);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x14000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 17);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00884022);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x01201048);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00040201);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 4), 4));
#endif

#ifdef CASE_T4V0
	//T4V0
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x0158737d);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 21);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x08000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 4);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000882);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x0000313c);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 13);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x07ffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 27));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));


	c = _mm_set1_epi32(0x00078c82);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x50000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 28);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x06000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 21);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x20a00000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 15);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000641);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 3);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 22), 4));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 10), 8));
#endif

#ifdef CASE_T4V1
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x00044151);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x10000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 21);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00500000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 10);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x50054000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 13);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x01000400);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 8);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000015);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 9);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x0007ffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 19));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));

	c = _mm_set1_epi32(0xef2bbcae);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x08000000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 27);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x2000a800);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00100200);
	t = _mm_and_si128(m[1], c);
	res[0] = _mm_or_si128(res[0], t);
	c = _mm_set1_epi32(0x00421100);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 6);
	res[0] = _mm_or_si128(res[0], s);

	c = _mm_set1_epi32(0x000000ca);
	res[1] = _mm_and_si128(m[1], c);
	c = _mm_set1_epi32(0x06a00000);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 21);
	res[1] = _mm_or_si128(res[1], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[1] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[1], _mm_bsrli_si128(res[1], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 6), 4));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 26), 8));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[1], 6), 8));


#endif

#ifdef CASE_T5V0
	//T5V0
	__m128i c, t, s, mAC, mAG, mMask, res[2];
	c = _mm_set1_epi32(0x00007d63);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x80000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 29);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000030);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00002080);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 2);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00001e01);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 7);
	res[0] = _mm_or_si128(res[0], s);

	mMask = _mm_set_epi32(0, 0, 0, 0x000fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	res[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 20));
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));


	c = _mm_set1_epi32(0x03ff829c);
	res[0] = _mm_and_si128(m[0], c);
	c = _mm_set1_epi32(0x20000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 29);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x5c000000);
	t = _mm_and_si128(m[0], c);
	s = _mm_srli_epi32(t, 16);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x00000044);
	t = _mm_and_si128(m[1], c);
	s = _mm_srli_epi32(t, 1);
	res[0] = _mm_or_si128(res[0], s);
	c = _mm_set1_epi32(0x0000010a);
	t = _mm_and_si128(m[1], c);
	s = _mm_slli_epi32(t, 5);
	res[0] = _mm_or_si128(res[0], s);

	res[0] = _mm_and_si128(_mm_set_epi32(0, 0, 0, 0xffffffff), _mm_or_si128(res[0], _mm_bsrli_si128(res[0], 8)));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_slli_epi32(res[0], 8), 4));
	vPattern[0] = _mm_or_si128(vPattern[0], _mm_bslli_si128(_mm_srli_epi32(res[0], 24), 8));
#endif

#ifdef CASE_C16
	//C16
	__m128i mAC, mAG, mMask;

	mMask = _mm_set_epi32(0, 0, 0, 0x0000ffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8)));
	vPattern[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 16));
#endif

#ifdef CASE_C20
	//C20
	__m128i mAC, mAG, mMask;

	mMask = _mm_set_epi32(0, 0, 0, 0x000fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8)));
	vPattern[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 20));
#endif

#ifdef CASE_C24
	//C24
	__m128i mAC, mAG, mMask;

	mMask = _mm_set_epi32(0, 0, 0, 0x00ffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8)));
	vPattern[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 24));
#endif

#ifdef CASE_C28
	//C28
	__m128i mAC, mAG, mMask;

	mMask = _mm_set_epi32(0, 0, 0, 0x0fffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8)));
	vPattern[0] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 28));
#endif

#ifdef CASE_C32
	//C32
	__m128i mAC, mAG, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	vPattern[0] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));
#endif

#ifdef CASE_C36
	//C36
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x0000000f);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 4));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C40
	//C40
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x000000ff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 8));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C44
	//C44
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x00000fff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 12));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C48
	//C48
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x0000ffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 16));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C52
	//C52
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x000fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 20));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C56
	//C56
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x00ffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 24));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

#ifdef CASE_C60
	//C60
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mMask = _mm_set_epi32(0, 0, 0, 0x0fffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8)));
	res[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 28));
	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif


#ifdef CASE_C64
	//C64
	__m128i mAC, mAG, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);
#endif

	return 0;
}


int patternReadSeed::comparePatternLibrary(char* vcPattern, char* vcRecord) {
	int ik;
	unsigned int uValue1, uValue2;

	for (int i = 1; i <= nSearchWhole; i++) {
		ik = nSearchBytes - 4 * i;
		uValue1 = *(unsigned int*)(vcPattern + ik);
		uValue2 = *(unsigned int*)(vcRecord + ik);
		if (uValue1 == uValue2) continue;
		if (uValue1 < uValue2) {
			return 1;
		}
		return -1;
	}
	uValue1 = *(unsigned int*)(vcPattern);
	uValue2 = *(unsigned int*)(vcRecord);
	if (uValue1 == uValue2) {
		return 0;
	}
	if (uValue1 < uValue2) {
		return 1;
	}

	return -1;
}

int patternReadSeed::findLibraryList(char* vc, int indFirst, int indLast, int* indListFirst, int* indListLast) {
	int iValueStart, iValueEnd, iValue3;
	int ind1, ind2, ind3;
	int ind1U, ind2L, ind1C, ind2C;
	char* vcPattern;
	bool isFound;
	*indListFirst = indFirst;
	*indListLast = indLast;

	if (indLast == indFirst) return 0;

	vcPattern = (char*)vPatternShift;

	iValueStart = comparePatternLibrary(vcPattern, vc + indFirst * recordSize);
	if (iValueStart > 0) {
		//printf("Start issue\n");
		*indListLast = indFirst;
		return 0;
	}

	iValueEnd = comparePatternLibrary(vcPattern, vc + (indLast - 1) * recordSize);
	if (iValueEnd < 0) {
		//printf("End issue\n");
		*indListFirst = indLast;
		return 0;
	}

	ind1 = indFirst;
	ind2 = indLast - 1;

	if (iValueStart == 0) {
		ind3 = indFirst;
	}
	else if (iValueEnd == 0) {
		ind3 = indLast - 1;
	}
	else {
		isFound = false;
		while (ind2 - ind1 > 1) {
			ind3 = (ind1 + ind2) / 2;
			iValue3 = comparePatternLibrary(vcPattern, vc + ind3 * recordSize);
			//printf("(%i, %i)\tpos: %i\tvalue :%i\n", ind1, ind2, ind3, iValue3);
			if (iValue3 == 0) {
				isFound = true;
				break;
			}
			if (iValue3 < 0) {
				ind1 = ind3;
			}
			else {
				ind2 = ind3;
			}
		}
		if (!isFound) {
			*indListLast = indFirst;
			return 0;
		}
		else {
			*indListFirst = ind3;
			*indListLast = ind3 + 1;
		}
	}

	if (ind3 == indFirst) {
		ind1 = ind3;
	}
	else {
		ind1U = ind3;
		while (ind1U - ind1 > 1) {
			ind1C = (ind1 + ind1U) / 2;
			iValue3 = comparePatternLibrary(vcPattern, vc + ind1C * recordSize);
			if (iValue3 == 0) {
				ind1U = ind1C;
			}
			else {
				ind1 = ind1C;
			}
		}
		ind1 = ind1U;
	}

	if (ind3 == indLast - 1) {
		ind2 = ind3;
	}
	else {
		ind2L = ind3;
		while (ind2 - ind2L > 1) {
			ind2C = (ind2 + ind2L) / 2;
			iValue3 = comparePatternLibrary(vcPattern, vc + ind2C * recordSize);
			if (iValue3 == 0) {
				ind2L = ind2C;
			}
			else {
				ind2 = ind2C;
			}
		}
		ind2 = ind2L;
	}

	//printf("Down: %i, up: %i\n", ind1, ind2);
	*indListFirst = ind1;
	*indListLast = ind2 + 1;

	return 0;
}

int patternReadSeed::formMaskRef() {
	unsigned int uu;
	int ik;
	for (int i = 0; i < 32; i++) {
		vOne[i] = 1 << i;
	}

	for (int k = 0; k < nLen; k++) {
		if (k % 32 == 0) uu = 0;
		if (vSeed[k] == '#' || vSeed[k] == '@' || vSeed[k] == '1' || vSeed[k] == '2') uu |= vOne[k % 32];
		if (k % 32 == 31 || k == nLen - 1) {
			ik = k / 32;
			vRef[ik] = _mm_set1_epi32(uu);
			vMask[ik] = _mm_set_epi32(0, 0, 0, uu);
		}
	}

	return 0;
}

int patternReadSeed::findWeight() {
	nDoubleWeight = 0;
	for (int i = 0; i < nLen; i++) {
		if (vSeed[i] == '1' || vSeed[i] == '#') nDoubleWeight += 2;
		if (vSeed[i] == '2' || vSeed[i] == '@') nDoubleWeight++;
	}

	recordSize = nDoubleWeight - 4 * nLetters + 32;
	if (recordSize % 8 > 0)recordSize += 8;
	recordSize /= 8;

	nSearchBits = nDoubleWeight - 4 * nLetters - nIndexLevel;
	nSearchBytes = nSearchBits / 8;
	if (nSearchBits % 8 > 0) nSearchBytes++;
	nSearchRem = nSearchBytes % 4;
	nSearchWhole = nSearchBytes / 4;
	if (nSearchRem == 0) {
		nSearchWhole--;
		nSearchRem = 4;
	}

	nPattern = nDoubleWeight / 128;
	if (nDoubleWeight % 128 > 0) nPattern++;
#ifdef WIN32
	vPattern = (__m128i*)_aligned_malloc(sizeof(__m128i) * nPattern, sizeof(__m128i));
	vPatternShift = (__m128i*)_aligned_malloc(sizeof(__m128i) * nPattern, sizeof(__m128i));
	vQ = (char*)_aligned_malloc(sizeof(__m128i) * nPattern, sizeof(__m128i));
#else
	vPattern = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nPattern);
	vPatternShift = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nPattern);
	vQ = (char*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nPattern);
#endif

	return 0;
}

patternReadSeed::patternReadSeed(char* vSeed_in, int nLetters_in, int nIndexLevel_in) {
	nLen = strlen(vSeed_in);
	strcpy(vSeed, vSeed_in);
	nLen128 = nLen / 32;
	if (nLen % 32 > 0) nLen128++;

	nLetters = nLetters_in;
	nIndexLevel = nIndexLevel_in;

	findWeight();

#ifdef WIN32
	vm = (__m128i*)_aligned_malloc(sizeof(__m128i) * nLen128, sizeof(__m128i));
	vMask = (__m128i*)_aligned_malloc(sizeof(__m128i) * nLen128, sizeof(__m128i));
	vRef = (__m128i*)_aligned_malloc(sizeof(__m128i) * nLen128, sizeof(__m128i));
	vOut = (__m128i*)_aligned_malloc(sizeof(__m128i) * (nLen128 + 1), sizeof(__m128i));
#else
	vm = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nLen128);
	vMask = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nLen128);
	vRef = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nLen128);
	vOut = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * (nLen128 + 1));

#endif

	mZero = _mm_set1_epi32(0);

	shiftBitFile = 4 * nLetters;

	uFileMask = 0;
	for (int i = 0; i < nLetters; i++) {
		uFileMask |= (15 << (4 * i));
	}

	formMaskRef();
}

int patternReadSeed::copyData(__m128i* v_in, int iShift_in) {
	int iRem, iWhole;
	iShift = iShift_in;
	iRem = iShift & 31;
	iWhole = iShift >> 5;

	m3 = mZero;

	for (int k = 0; k < nLen128; k++) {
		m1 = _mm_srli_epi32(v_in[iWhole + k], iRem);
		m2 = _mm_slli_epi32(v_in[iWhole + k + 1], 32 - iRem);
		vm[k] = _mm_and_si128(vRef[k], _mm_or_si128(m1, m2));

		m4 = _mm_or_si128(vm[k], _mm_bsrli_si128(vm[k], 8));
		m5 = _mm_or_si128(m4, _mm_bsrli_si128(m4, 4));
		m3 = _mm_or_si128(m3, _mm_xor_si128(m5, vRef[k]));
	}

	isOk = (_mm_extract_epi32(m3, 0) == 0);
	if (isOk) {
		formSignature();
		return 0;
	}
	return -1;
}

int patternReadSeed::formSignature() {
	formSignatureM(vm, vPattern);
	findFix();

	return 0;
}

int patternReadSeed::findFix() {
	__m128i m1, m2, m3, m4;
	int n128;
	int iValue, i1, i2, ii1, ii2, ij1, ij2, ir;

	indFile = uFileMask & _mm_extract_epi32(vPattern[0], 0);

	for (int i = 0; i < nPattern; i++) {
		m1 = _mm_srli_epi64(vPattern[i], shiftBitFile);
		m2 = _mm_slli_epi64(_mm_bsrli_si128(vPattern[i], 8), 64 - shiftBitFile);
		vPatternShift[i] = _mm_or_si128(m1, m2);
	}

	for (int i = 1; i < nPattern; i++) {
		m2 = _mm_slli_epi64(_mm_bslli_si128(vPattern[i], 8), 64 - shiftBitFile);
		vPatternShift[i - 1] = _mm_or_si128(vPatternShift[i - 1], m2);
	}


	i2 = nDoubleWeight - 1;
	i1 = i2 - nIndexLevel + 1;

	ir = i1 % 32;
	ii1 = i1 / 128;
	ij1 = (i1 % 128) / 32;
	ii2 = i2 / 128;
	ij2 = (i2 % 128) / 32;

	if (ij1 == 0) m1 = _mm_bsrli_si128(vPattern[ii1], 0);
	if (ij1 == 1) m1 = _mm_bsrli_si128(vPattern[ii1], 4);
	if (ij1 == 2) m1 = _mm_bsrli_si128(vPattern[ii1], 8);
	if (ij1 == 3) m1 = _mm_bsrli_si128(vPattern[ii1], 12);

	if (ij2 == 0) m2 = _mm_srli_si128(vPattern[ii2], 0);
	if (ij2 == 1) m2 = _mm_srli_si128(vPattern[ii2], 4);
	if (ij2 == 2) m2 = _mm_srli_si128(vPattern[ii2], 8);
	if (ij2 == 3) m2 = _mm_srli_si128(vPattern[ii2], 12);

	m3 = _mm_unpacklo_epi32(m1, m2);
	m4 = _mm_srli_epi64(_mm_slli_epi64(m3, 64 - ir - nIndexLevel), 64 - nIndexLevel);
	indIndex = _mm_extract_epi32(m4, 0);


	iValue = nDoubleWeight - shiftBitFile;

	n128 = iValue / 128;
	iValue -= n128 * 128;

	return 0;
}

patternReadSeed::~patternReadSeed() {
#ifdef WIN32
	if (vm != nullptr) { _aligned_free(vm); vm = nullptr; }
	if (vRef != nullptr) { _aligned_free(vRef); vRef = nullptr; }
	if (vMask != nullptr) { _aligned_free(vMask); vMask = nullptr; }
	if (vPattern != nullptr) { _aligned_free(vPattern); vPattern = nullptr; }
	if (vPatternShift != nullptr) { _aligned_free(vPatternShift); vPatternShift = nullptr; }
	if (vOut != nullptr) { _aligned_free(vOut); vOut = nullptr; }
	if (vQ != nullptr) { _aligned_free(vQ); vQ = nullptr; }
#else
	if (vm != nullptr) { free(vm); vm = nullptr; }
	if (vRef != nullptr) { free(vRef); vRef = nullptr; }
	if (vMask != nullptr) { free(vMask); vMask = nullptr; }
	if (vPattern != nullptr) { free(vPattern); vPattern = nullptr; }
	if (vPatternShift != nullptr) { free(vPatternShift); vPatternShift = nullptr; }
	if (vOut != nullptr) { free(vOut); vOut = nullptr; }
	if (vQ != nullptr) { free(vQ); vQ = nullptr; }
#endif
}

contP::contP() {
	voBest = nullptr;
	voCountCandidates = nullptr;
	voCountSolutions = nullptr;
	nThreads = 0;
	nPatterns = 0;
	nLibraryEntries = 0;


	pBD = nullptr;


	uChromStart = nullptr;

	pvPositions = nullptr;
	vCountPositions = nullptr;
	vnLibraryRecords = nullptr;
	pvLibraryData = nullptr;
	pvLibraryStart = nullptr;

	vmRef = nullptr;
	fiLeft = nullptr;
	fiRight = nullptr;
	fiACGT = nullptr;
	fo = nullptr;
	foStat = nullptr;
	fi = nullptr;
}

contP::~contP() {
	if (pBD != nullptr) {
		for (int i = 0; i < nThreads; i++) {
			delete (pBD[i]); pBD[i] = nullptr;
			//if (pBD[i] != nullptr) { delete (pBD[i]); pBD[i] = nullptr; }
		}
		free(pBD); pBD = nullptr;
	}
	if (voBest != nullptr) { free(voBest); voBest = nullptr; }
	if (voCountCandidates != nullptr) { free(voCountCandidates); voCountCandidates = nullptr; }
	if (voCountSolutions != nullptr) { free(voCountSolutions); voCountSolutions = nullptr; }

	if (uChromStart != nullptr) { free(uChromStart); uChromStart = nullptr; }

	if (vnLibraryRecords != nullptr) {
		free(vnLibraryRecords); vnLibraryRecords = nullptr;
	}
	if (pvLibraryData != nullptr) {
		for (int i = 0; i < nLibraryFiles; i++) {
			if (pvLibraryData[i] != nullptr) { free(pvLibraryData[i]); pvLibraryData[i] = nullptr; }
		}
		free(pvLibraryData); pvLibraryData = nullptr;
	}
	if (pvLibraryStart != nullptr) {
		for (int i = 0; i < nLibraryFiles; i++) {
			if (pvLibraryStart[i] != nullptr) { free(pvLibraryStart[i]); pvLibraryStart[i] = nullptr; }
		}
		free(pvLibraryStart); pvLibraryStart = nullptr;
	}

	if (pvPositions != nullptr) {
		for (int i = 0; i < nPatterns; i++) {
			if (pvPositions[i] != nullptr) {
				free(pvPositions[i]); pvPositions[i] = nullptr;
			}
		}
		free(pvPositions); pvPositions = nullptr;
	}
	if (vCountPositions != nullptr) {
		free(vCountPositions); vCountPositions = nullptr;
	}

#ifdef WIN32
	if (vmRef != nullptr) { _aligned_free(vmRef); vmRef = nullptr; }
#else
	if (vmRef != nullptr) { free(vmRef); vmRef = nullptr; }
#endif
	if (fiLeft != nullptr) { fclose(fiLeft); fiLeft = nullptr; }
	if (fiRight != nullptr) { fclose(fiRight); fiRight = nullptr; }
	if (fiACGT != nullptr) { fclose(fiACGT); fiACGT = nullptr; }
	if (fo != nullptr) { fclose(fo); fo = nullptr; }
	if (foStat != nullptr) { fclose(foStat); foStat = nullptr; }
	if (fi != nullptr) { fclose(fi); fi = nullptr; }
}

int blockData::measureSimilarity(unsigned int upos, __m128i* vR, int* countT, int* countV) {
	__m128i mRef, * vRLoc, m1, m2, m3, mShuf;
	int refBlock, refRem;
	refBlock = upos >> 5;
	refRem = upos & 31;
	vRLoc = vmRef + refBlock;
	*countT = 0;
	*countV = 0;
	for (int i = 0; i < readMsize; i++) {
		mRef = _mm_or_si128(_mm_srli_epi32(vRLoc[i], refRem), _mm_slli_epi32(vRLoc[i + 1], 32 - refRem));
		m1 = _mm_and_si128(mSim[i], _mm_and_si128(mRef, vR[i]));
		m2 = _mm_or_si128(m1, _mm_bsrli_si128(m1, 8));
		m3 = _mm_xor_si128(mSim[i], _mm_or_si128(m2, _mm_bsrli_si128(m2, 4)));
		* countV += _popcnt32(_mm_extract_epi32(m3, 0));

		mShuf = _mm_shuffle_epi32(mRef, 78);
		m1 = _mm_and_si128(mSim[i], _mm_and_si128(mShuf, vR[i]));
		m2 = _mm_or_si128(m1, _mm_bsrli_si128(m1, 8));
		m3 = _mm_or_si128(m2, _mm_bsrli_si128(m2, 4));
		* countT += _popcnt32(_mm_extract_epi32(m3, 0));
	}

	*countV -= *countT;

	return 0;
}

int contP::getInfoReferenceData() {
	char mStr[1000];
	fi = fopen(inputInfoReference, "r");
	if (fi == nullptr) {
		printf("Error: cannot open file \"%s\"\n", inputInfoReference);
		return -1;
	}

	nChromosomes = 0;
	while (fscanf(fi, "%[^\n]\n", mStr) > 0) {
		nChromosomes++;
	}
	nChromosomes /= 3;
	printf("Number of blocks in chrom data: %i\n", nChromosomes);

	uChromStart = (long long*)malloc(sizeof(long long) * nChromosomes);

	rewind(fi);
	nChromosomes = 0;
	while (fscanf(fi, "%[^\n]\n", mStr) > 0) {
		if (nChromosomes % 3 == 1) {
			uChromStart[nChromosomes / 3] = (long long)2 * atoll(mStr);
			printf("%llu\n", uChromStart[nChromosomes / 3]);
		}
		nChromosomes++;
	}
	nChromosomes /= 3;
	fclose(fi); fi = nullptr;
	return 0;
}

int blockData::formStringPosition(unsigned int upos, char* vC) {
	long long lpos;
	lpos = (long long)upos;

	for (int i = nChromosomes - 1; i >= 0; i--) {
		if (lpos >= uChromStart[i]) {
			sprintf(vC, "(%i)%lli", i + 1, (long long)(1) + lpos - uChromStart[i]);
			break;
		}
	}

	return 0;
}

int blockData::findPairedList() {
	unsigned int* v1, * v2, * vout;
	unsigned int minRange, maxRange;
	int n1, n2, nout;
	int shiftMin, shiftMax;


	for (int ik = 0; ik < 2; ik++) {
		nout = 0;
		if (ik == 0) {
			v1 = pvMergedCandidatePosition[0];
			v2 = pvMergedCandidatePosition[3];
			shiftMin = minPairDistance;
			shiftMax = maxPairDistance;
			n1 = vCountMergedCandidate[0];
			n2 = vCountMergedCandidate[3];
			vout = vPairPosition03;
		}
		else {
			v1 = pvMergedCandidatePosition[1];
			v2 = pvMergedCandidatePosition[2];
			shiftMin = minPairDistance;
			shiftMax = maxPairDistance;
			n1 = vCountMergedCandidate[1];
			n2 = vCountMergedCandidate[2];
			vout = vPairPosition12;
		}

		for (int i = 0; i < n1; i++) {
			minRange = v1[i] + shiftMin;
			maxRange = v1[i] + shiftMax;
			for (int j = 0; j < n2; j++) {
				if (v2[j] > maxRange) break;
				if (v2[j] < minRange) continue;
				vout[2 * nout] = v1[i];
				vout[2 * nout + 1] = v2[j];
				nout++;
			}
		}
		if (ik == 0) {
			nCountPair03 = nout;
		}
		else {
			nCountPair12 = nout;
		}
	}

	return 0;
}

int blockData::mergeCandiateListSpaced(int ind) {
	int countLocal, nTotal, sumTotal;
	unsigned int valueMin, valueMax;

	for (int i = 0; i < nPatterns; i++) {
		vIndexCurrent[i] = 0;
	}

	valueMin = 0xffffffff;
	valueMax = 0;

	sumTotal = 0;

	for (int i = 0; i < nPatterns; i++) {
		countLocal = vCountCandidate[nPatterns * ind + i];
		sumTotal += countLocal;
		if (countLocal == 0) continue;
		if (valueMin > pvCandidatePosition[nPatterns * ind + i][0]) valueMin = pvCandidatePosition[nPatterns * ind + i][0];
		if (valueMax < pvCandidatePosition[nPatterns * ind + i][countLocal - 1]) valueMax = pvCandidatePosition[nPatterns * ind + i][countLocal - 1];
	}

	if (sumTotal == 0) {
		vCountMergedCandidate[ind] = 0;
		return 0;
	}

	pvMergedCandidatePosition[ind][0] = valueMin;
	nTotal = 1;


	if (valueMin < valueMax) {
		do {
			for (int i = 0; i < nPatterns; i++) {
				countLocal = vCountCandidate[nPatterns * ind + i];
				if (vIndexCurrent[i] >= countLocal)continue;
				if (pvCandidatePosition[nPatterns * ind + i][vIndexCurrent[i]] == valueMin) {
					vIndexCurrent[i]++;
				}
			}

			valueMin = valueMax;
			for (int i = 0; i < nPatterns; i++) {
				countLocal = vCountCandidate[nPatterns * ind + i];
				if (vIndexCurrent[i] >= countLocal)continue;
				if (pvCandidatePosition[nPatterns * ind + i][vIndexCurrent[i]] < valueMin) {
					valueMin = pvCandidatePosition[nPatterns * ind + i][vIndexCurrent[i]];
				}
			}
			pvMergedCandidatePosition[ind][nTotal] = valueMin;
			nTotal++;
			if (valueMax == valueMin) break;
		} while (true);
	}

	vCountMergedCandidate[ind] = nTotal;

	return 0;
}

int contP::printLibraryRecord(char* vc, int ind) {
	char* vcLoc;
	vcLoc = vc + ind * recordSize;
	fprintf(fo, "%8i\t", ind);
	for (int j = 0; j < recordSize - 4; j++) {
		fprintf(fo, "%1x%1x", vcLoc[j] & 15, 15 & (vcLoc[j] >> 4));
	}
	fprintf(fo, "\t%012u\n", *(unsigned int*)(vcLoc + recordSize - 4));
	return 0;
}

int contP::loadLibrary() {
	recordSize = seedDoubleWeight + 32 - 4 * seedNLetters;
	if (recordSize % 8 > 0) recordSize += 8;
	recordSize /= 8;

	printf("Size of records: %i\n", recordSize);

	nLibraryFiles = 1 << (4 * seedNLetters);
	printf("Number of library files: %i\n", nLibraryFiles);

	vnLibraryRecords = (int*)malloc(sizeof(int) * nLibraryFiles);

	pvLibraryData = (char**)malloc(sizeof(char*) * nLibraryFiles);
	for (int i = 0; i < nLibraryFiles; i++) {
		pvLibraryData[i] = nullptr;
	}
	pvLibraryStart = (int**)malloc(sizeof(int*) * nLibraryFiles);
	for (int i = 0; i < nLibraryFiles; i++) {
		pvLibraryStart[i] = nullptr;
	}

	sprintf(libraryTemplate, "%s/sorted/%%0%ix.bin", inputLibrary, seedNLetters);
	printf("Library template: %s\n", libraryTemplate);

	nLibraryEntries = 1 << seedIndexLevel;

	printf("Library entries: %i\n", nLibraryEntries);


	for (int i = 0; i < nLibraryFiles; i++) {
		pvLibraryStart[i] = (int*)malloc(sizeof(int) * (nLibraryEntries + 1));
		pvLibraryStart[i][0] = 0;
	}

	for (int i = 0; i < nLibraryFiles; i++) {
		pvLibraryData[i] = nullptr;
	}

	int *vLS;
	int nLR;
	char* vC;

	for (int k = 0; k < nLibraryFiles; k++) {
		if (k % 1000 == 0) printf("Library files: %i\n", k);
		sprintf(inputFile, libraryTemplate, k);
		fi = fopen(inputFile, "rb");
		if (fi == nullptr) {
			printf("Error: cannot open library file %s\n", inputFile);
			return -1;
		}

		vLS = pvLibraryStart[k];
		fread(vLS + 1, sizeof(int), nLibraryEntries, fi);
		nLR = vLS[nLibraryEntries];
		vC = (char*)malloc(sizeof(char) * recordSize * nLR);
		fread(vC, sizeof(char), recordSize * nLR, fi);
		fclose(fi); fi = nullptr;
		pvLibraryData[k] = vC;
		vnLibraryRecords[k] = nLR;
	}

	return 0;
}

int contP::printInfo(){
	printf("Usage\n");
	printf("1) FQB file (first)\n");
	printf("2) FQB file (second)\n");
	printf("3) path to library folder\n");
	printf("4) reference acgt file\n");
	printf("5) reference iacgt file\n");
	printf("6) output folder\n");
	printf("7) distance (min)\n");
	printf("8) distance (max)\n\n");

	return 0;
}

int contP::setParameters() {
	int iPrint, iRem;
	iPrint = 10;
	iRem = PRINT_READ_SHIFT - iPrint - 4;
	sprintf(printReadTemplateStart, "%%%ii", iPrint);

	for (int i = 0; i < iRem; i++) {
		printReadTemplateGap[i] = ' ';
	}
	printReadTemplateGap[iRem] = '\0';

	iPrint = 10;
	iRem = PRINT_READ_SHIFT - iPrint - PRINT_EXTRA_REF;
	sprintf(printRefTemplateStart, "%%%ii", iPrint);

	for (int i = 0; i < iRem; i++) {
		printRefTemplateGap[i] = ' ';
	}
	printRefTemplateGap[iRem] = '\0';

	return 0;
}

int contP::printReference(unsigned int uPos) {
	//int nRLenPrint;
	//int nmStart, nmEnd;
	unsigned int u1, u2;
	int ind1, ind2, du, iRes, nSame;
	__m128i* vuLoc;
	__m128i m1, m2, m3, m4, mW, mMask;
	char* vStringLoc;

	mMask = _mm_set1_epi32(1);
	mW = _mm_set_epi32(8, 4, 2, 1);

	nLenRef = readLengthLeft + 2 * PRINT_EXTRA_REF;
	u1 = uPos - PRINT_EXTRA_REF;
	u2 = uPos + PRINT_EXTRA_REF + readLengthLeft - 1;
	ind1 = u1 / 32;
	ind2 = u2 / 32;
	du = ind2 - ind1 + 1;

	vuLoc = vmRef + ind1;

	vStringLoc = stringRef + u1 % 32;

	for (int k = 0; k < du; k++) {
		for (int j = 0; j < 32; j++) {
			m1 = _mm_and_si128(_mm_srli_epi32(vuLoc[0], j), mMask);
			m2 = _mm_madd_epi16(m1, mW);
			m3 = _mm_hadd_epi32(m2, m2);
			m4 = _mm_hadd_epi32(m3, m3);
			iRes = _mm_extract_epi32(m4, 0);
			switch (iRes) {
			case 1: stringRef[k * 32 + j] = 'a'; break;
			case 2: stringRef[k * 32 + j] = 'c'; break;
			case 4: stringRef[k * 32 + j] = 'g'; break;
			case 8: stringRef[k * 32 + j] = 't'; break;
			default: stringRef[k * 32 + j] = 'n'; break;
			}
		}
	}
	vStringLoc[nLenRef] = '\0';

	nSame = 0;

	for (int i = 0; i < readLengthLeft; i++) {
		if (vStringLoc[i + PRINT_EXTRA_REF] != stringRead[i]) {
			vStringLoc[i + PRINT_EXTRA_REF] -= 32;
		}
		else {
			nSame++;
		}
	}

	fprintf(fo, printRefTemplateStart, uPos);
	fprintf(fo, "%s", printRefTemplateGap);
	fprintf(fo, "%s%8i\n", vStringLoc, nSame);

	return 0;
}

int blockData::processReads(int indRelativeStart_in, int nReadsToProcess_in) {
	int ind, sumM, totalPairs;
	int indFirst, indLast, indListFirst, indListLast, indexF, posStart;
	int countT1, countT2, countV1, countV2;
	int iTV, iBestTV;
	long long indAbsolute;
	char cChrOut1[100], cChrOut2[100];
	char* vcLoc;

	long long ib1, ib2;

	unsigned int* uCP, uLoc;

	indRelativeStart = indRelativeStart_in;
	nReadsToProcess = nReadsToProcess_in;

	for (int ik = 0; ik < nReadsToProcess; ik++) {
		ind = indRelativeStart + ik;
		indAbsolute = (long long)ind + indReadStart + 1;
		fprintf(fo, ">%lli\n", indAbsolute);

		ib1 = (long long)(ind) * (long long)(blockMsize);
		ib2 = ib1 + (long long) (readMsize);

		for (int iq = 0; iq < 4; iq++) {
			switch (iq) {
				case 0: vmReadLoc = vmReadLeft + ib1; break;
				case 1: vmReadLoc = vmReadRight + ib1; break;
				case 2: vmReadLoc = vmReadLeft + ib2; break;
				case 3: vmReadLoc = vmReadRight + ib2; break;
			}

			for (int k = 0; k < nPatterns; k++) {
				uCP = pvCandidatePosition[iq * nPatterns + k];
				posStart = indPatternStartPosition[k];
				vCountCandidate[iq * nPatterns + k] = 0;
				if (pPRS[k]->copyData(vmReadLoc, posStart) == -1)continue;
				indexF = pPRS[k]->indFile;
				vcLoc = pvLibraryData[indexF];

				if (vcLoc == nullptr) continue;

				indFirst = pvLibraryStart[indexF][pPRS[k]->indIndex];
				indLast = pvLibraryStart[indexF][pPRS[k]->indIndex + 1];

				pPRS[k]->findLibraryList(vcLoc, indFirst, indLast, &indListFirst, &indListLast);

				vCountCandidate[iq * nPatterns + k] = indListLast - indListFirst;

				for (int i = indListFirst; i < indListLast; i++) {
					uLoc = *(unsigned int*)(vcLoc + (i + 1) * recordSize - 4);
					uCP[i - indListFirst] = (unsigned int)((long long)(uLoc) - (long long)(posStart));
				}
			}
			mergeCandiateListSpaced(iq);
		}

		totalPairs = 0;
		sumM = 0;
		for (int i = 0; i < 4; i++) {
			sumM += vCountMergedCandidate[i];
		}

		fprintf(fo, "MList\t%i\t%i\t%i\t%i\t%i\n", vCountMergedCandidate[0], vCountMergedCandidate[1], vCountMergedCandidate[2], vCountMergedCandidate[3], sumM);

		findPairedList();
		iBestTV = 999;

		totalPairs = nCountPair03 + nCountPair12;

		for (int i = 0; i < nCountPair03; i++) {
			formStringPosition(vPairPosition03[2 * i], cChrOut1);
			formStringPosition(vPairPosition03[2 * i + 1], cChrOut2);
			measureSimilarity(vPairPosition03[2 * i], vmReadLeft + ib1, &countT1, &countV1);
			measureSimilarity(vPairPosition03[2 * i + 1], vmReadRight + ib2, &countT2, &countV2);
			//fprintf(fo, "1\t%u\t%s\tV%i\tT%i\t", vPairPosition03[2 * i], cChrOut1, countV1, countT1);
			//fprintf(fo, "%u\t%s\tV%i\tT%i\t%i\n", vPairPosition03[2 * i + 1], cChrOut2, countV2, countT2, int((long long)vPairPosition03[2 * i + 1] - (long long)vPairPosition03[2 * i]));

			fprintf(fo, "1\t%s\tV%i\tT%i\t", cChrOut1, countV1, countT1);
			fprintf(fo, "%s\tV%i\tT%i\t%i\n", cChrOut2, countV2, countT2, int((long long)vPairPosition03[2 * i + 1] - (long long)vPairPosition03[2 * i]));

			iTV = 2 * (countV1 + countV2) + countT1 + countT2;
			if (iTV < iBestTV) iBestTV = iTV;
		}

		for (int i = 0; i < nCountPair12; i++) {
			formStringPosition(vPairPosition12[2 * i], cChrOut1);
			formStringPosition(vPairPosition12[2 * i + 1], cChrOut2);
			measureSimilarity(vPairPosition12[2 * i], vmReadRight + ib1, &countT1, &countV1);
			measureSimilarity(vPairPosition12[2 * i + 1], vmReadLeft + ib2, &countT2, &countV2);
			//fprintf(fo, "2\t%u\t%s\tV%i\tT%i\t", vPairPosition12[2 * i], cChrOut1, countV1, countT1);
			//fprintf(fo, "%u\t%s\tV%i\tT%i\t%i\n", vPairPosition12[2 * i + 1], cChrOut2, countV2, countT2, int((long long)vPairPosition12[2 * i + 1] - (long long)vPairPosition12[2 * i]));
			fprintf(fo, "2\t%s\tV%i\tT%i\t", cChrOut1, countV1, countT1);
			fprintf(fo, "%s\tV%i\tT%i\t%i\n", cChrOut2, countV2, countT2, int((long long)vPairPosition12[2 * i + 1] - (long long)vPairPosition12[2 * i]));
			iTV = 2 * (countV1 + countV2) + countT1 + countT2;
			if (iTV < iBestTV) iBestTV = iTV;
		}

		fprintf(fo, "<Best\t%i\n\n", iBestTV);

		voBestLoc[ind] = iBestTV;
		voCountCandidatesLoc[ind] = sumM;
		voCountSolutionsLoc[ind] = totalPairs;
	}
	return 0;
}

int contP::getSeedInfo() {
	sprintf(inputFile, "%s/info.txt", inputLibrary);
	printf("Library info file: \"%s\"\n", inputFile);
	fi = fopen(inputFile, "r");
	if (fi == nullptr) {
		printf("Error: cannot open the file \"%s\".\n", inputFile);
		return -1;
	}
	fscanf(fi, "%s\t%i\n", sTemp, &seedLength);
	printf("Length (seed): %i\n", seedLength);
	fscanf(fi, "%s\t%i\n", sTemp, &seedDoubleWeight);
	fscanf(fi, "%s\t%s\n", sTemp, sSeed);
	fscanf(fi, "%s\t%i\n", sTemp, &seedNLetters);
	fscanf(fi, "%s\t%i\n", sTemp, &seedIndexLevel);

	fclose(fi); fi = nullptr;

	printf("Here 2\n");

	printf("Length (seed): %i\n", seedLength);
	printf("Double weight: %i\n", seedDoubleWeight);
	printf("Seed: %s\n", sSeed);
	printf("Number of letters: %i\n", seedNLetters);
	printf("Index level: %i\n\n", seedIndexLevel);
	return 0;
}

int contP::getReadInfo() {
	//Left
	printf("Read file (left): %s\n", inputPairLeft);
	fiLeft = fopen(inputPairLeft, "rb");
	if (fiLeft == nullptr) {
		printf("Error: cannot open the file \"%s\".\n", inputPairLeft);
		return -1;
	}
	fread(&readLengthLeft, sizeof(int), 1, fiLeft);
	fread(&countReadLeft, sizeof(long long), 1, fiLeft);

#ifdef WIN32
	_fseeki64(fiLeft, 0, SEEK_END);
	fileSizeLeft = _ftelli64(fiLeft);
#else
	fseeko64(fiLeft, 0, SEEK_END);
	fileSizeLeft = ftello64(fiLeft);
#endif

	fclose(fiLeft); fiLeft = nullptr;

	printf("Read length (left): %i\n", readLengthLeft);
	printf("# reads (left): %lli\n", countReadLeft);
	printf("File size (left): %lli\n\n", fileSizeLeft);

	//Right
	printf("Read file (right): %s\n", inputPairRight);
	fiRight = fopen(inputPairRight, "rb");
	if (fiRight == nullptr) {
		printf("Error: cannot open the file \"%s\".\n", inputPairRight);
		return -2;
	}
	fread(&readLengthRight, sizeof(int), 1, fiRight);
	fread(&countReadRight, sizeof(long long), 1, fiRight);

#ifdef WIN32
	_fseeki64(fiRight, 0, SEEK_END);
	fileSizeRight = _ftelli64(fiRight);
#else
	fseeko64(fiRight, 0, SEEK_END);
	fileSizeRight = ftello64(fiRight);
#endif

	fclose(fiRight); fiRight = nullptr;

	printf("Read length (right): %i\n", readLengthRight);
	printf("# reads (right): %lli\n", countReadRight);
	printf("File size (right): %lli\n\n", fileSizeRight);

	if (fileSizeLeft != fileSizeRight) {
		printf("Error: files have different sizes.\n");
		return -3;
	}

	if (readLengthLeft != readLengthRight) {
		printf("Error: lengths of reads are not the same for two files\n");
		return -4;
	}

	if (countReadLeft != countReadRight) {
		printf("Error: number of reads are different for two files.\n");
		return -5;
	}

	readBinarySize = readLengthLeft;
	if (readLengthLeft % 32 > 0) readBinarySize += 32;
	readBinarySize /= 32;
	readBinarySize *= 16;

	printf("Size of data for each read: %i\n", readBinarySize);

	readBlockSize = 2 * readBinarySize;

	if ((long long)(readBlockSize)*countReadLeft + 12 != fileSizeLeft) {
		printf("Error: wrong number of reads.\n");
		return -6;
	}

	readMsize = readBinarySize / 16;
	blockMsize = readBlockSize / 16;

	return 0;
}


int contP::getReferenceData() {
	long long fs;

	printf("Reference file: %s\n", inputReference);
	fiACGT = fopen(inputReference, "rb");
	if (fiACGT == nullptr) {
		printf("Error: cannot open reference file \"%s\"\n", inputReference);
		return -1;
	}

#ifdef WIN32
	_fseeki64(fiACGT, 0, SEEK_END);
	fs = _ftelli64(fiACGT);
#else
	fseeko64(fiACGT, 0, SEEK_END);
	fs = ftello64(fiACGT);
#endif

	printf("File size (reference): %lli\n", fs);
	nLenRef = (int)(fs / 16);

	rewind(fiACGT);

#ifdef WIN32
	vmRef = (__m128i*)_aligned_malloc(sizeof(__m128i) * nLenRef, sizeof(__m128i));
#else
	vmRef = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * nLenRef);
#endif
	if (vmRef == nullptr) {
		printf("Error: cannot allocate memory\n");
		return -2;
	}

	fread(vmRef, sizeof(__m128i), nLenRef, fiACGT);

	fclose(fiACGT); fiACGT = nullptr;
	printf("Reference data: done\n\n");

	return 0;
}

int blockData::formPatternStartPositions() {
	indPatternStartPosition = (int*)malloc(sizeof(int) * nPatterns);
	for (int i = 0; i < nPatterns; i++) {
		indPatternStartPosition[i] = i;
	}

	vCountCandidate = (int*)malloc(sizeof(int) * 4 * nPatterns);
	pvCandidatePosition = (unsigned int**)malloc(sizeof(unsigned int*) * 4 * nPatterns);
	for (int i = 0; i < 4 * nPatterns; i++) {
		void *p;
		p = malloc(sizeof(unsigned int) * MAX_RECORDS);
		pvCandidatePosition[i] = (unsigned int*)p;
	}

	pvMergedCandidatePosition = (unsigned int**)malloc(sizeof(unsigned int*) * 4);
	for (int i = 0; i < 4; i++) {
		pvMergedCandidatePosition[i] = (unsigned int*)malloc(sizeof(unsigned int) * MAX_RECORDS);
	}

	vIndexCurrent = (int*)malloc(sizeof(int) * nPatterns);

	vPairPosition03 = (unsigned int*)malloc(sizeof(unsigned int) * 2 * MAX_RECORDS);
	vPairPosition12 = (unsigned int*)malloc(sizeof(unsigned int) * 2 * MAX_RECORDS);

	return 0;
}

int contP::startProcessing(int nArg, char **vArg) {
	int nReads2read, nReads, nW, nR, uQ;
	long long rStart;
	long long posLL, posRR;

	if(nArg != 9){
		printf("Error: wrong number of arguments\n\n");
		printInfo();
		return -1;
	}

	sprintf(inputPairLeft, "%s", vArg[1]);
	sprintf(inputPairRight, "%s", vArg[2]);
	sprintf(inputLibrary, "%s", vArg[3]);
	sprintf(inputReference, "%s", vArg[4]);
	sprintf(inputInfoReference, "%s", vArg[5]);
	sprintf(outputFolder, "%s", vArg[6]);
	minPairDistance = atoi(vArg[7]);
	maxPairDistance = atoi(vArg[8]);

	if (setParameters() != 0) return -1;
	if (getSeedInfo() != 0) return -2;
	if (getReadInfo() != 0) return -3;
	if (getInfoReferenceData() != 0) return -4;
	if (getReferenceData() != 0) return -5;
	if (loadLibrary() != 0) return -6;

	nPatterns = readLengthLeft - seedLength + 1;

	fiLeft = fopen(inputPairLeft, "rb");
	fiRight = fopen(inputPairRight, "rb");

	nReads = countReadLeft;
	nReads2read = MAX_READ_COUNT;

	voBest = (int*)malloc(sizeof(int) * nReads);
	voCountCandidates = (int*)malloc(sizeof(int) * nReads);
	voCountSolutions = (int*)malloc(sizeof(int) * nReads);

	nW = nReads / nReads2read;
	nR = nReads % nReads2read;

#pragma omp parallel
	{
		nThreads = omp_get_num_threads();
	}

	pBD = (blockData**)malloc(sizeof(blockData*) * nThreads);

	for (int i = 0; i < nThreads; i++) {
		pBD[i] = new blockData();

#ifdef WIN32
		pBD[i]->vmReadLeft = (__m128i*)_aligned_malloc(sizeof(__m128i) * (blockMsize * MAX_READ_COUNT + 1), sizeof(__m128i));
		pBD[i]->vmReadRight = (__m128i*)_aligned_malloc(sizeof(__m128i) * (blockMsize * MAX_READ_COUNT + 1), sizeof(__m128i));
#else
		pBD[i]->vmReadLeft = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * (blockMsize * MAX_READ_COUNT + 1));
		pBD[i]->vmReadRight = (__m128i*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * (blockMsize * MAX_READ_COUNT + 1));
#endif
		pBD[i]->blockMsize = blockMsize;
		pBD[i]->readMsize = readMsize;
		pBD[i]->voBestLoc = (int *)malloc(sizeof(int)*MAX_READ_COUNT);
		pBD[i]->voCountCandidatesLoc = (int *)malloc(sizeof(int)*MAX_READ_COUNT);
		pBD[i]->voCountSolutionsLoc = (int *)malloc(sizeof(int)*MAX_READ_COUNT);
		pBD[i]->nPatterns = nPatterns;
		pBD[i]->pvLibraryData = pvLibraryData;
		pBD[i]->pvLibraryStart = pvLibraryStart;
		pBD[i]->vmRef = vmRef;
		pBD[i]->minPairDistance = minPairDistance;
		pBD[i]->maxPairDistance = maxPairDistance;
		pBD[i]->recordSize = recordSize;
		pBD[i]->nChromosomes = nChromosomes;

		pBD[i]->uChromStart = (long long*)malloc(sizeof(long long) * nChromosomes);
		for (int ii = 0; ii < nChromosomes; ii++) {
			pBD[i]->uChromStart[ii] = uChromStart[ii];
		}

		pBD[i]->pPRS = (patternReadSeed**)malloc(sizeof(patternReadSeed*) * nPatterns);
		for (int ii = 0; ii < nPatterns; ii++) {
			pBD[i]->pPRS[ii] = new patternReadSeed(sSeed, seedNLetters, seedIndexLevel);
		}

#ifdef WIN32
		pBD[i]->mSim = (__m128i*) _aligned_malloc(sizeof(__m128i) * readMsize, sizeof(__m128i));
#else
		pBD[i]->mSim = (__m128i*) aligned_alloc(sizeof(__m128i), sizeof(__m128i) * readMsize);
#endif
		for (int ii = 0; ii < readMsize; ii++) {
			uQ = 0;
			for (int j = 0; j < 32; j++) {
				if (32 * ii + j < readLengthLeft) uQ |= (1 << j);
			}
			pBD[i]->mSim[ii] = _mm_set1_epi32(uQ);
		}

		pBD[i]->formPatternStartPositions();
	}

	int *vState;
	vState = (int *)malloc(sizeof(int)*nW);
	for(int i = 0; i < nW; i++){
		vState[i] = 0;
	}

#pragma omp parallel
	{
		blockData* bdLoc;
		int tid, ib, iFO;
		long long rStartLoc, posL, posR;
		char outputFileLoc[1000];
		tid = omp_get_thread_num();
		bdLoc = pBD[tid];

//#pragma omp for
		for (ib = 0; ib < nW; ib++) {
			iFO = 0;
#pragma omp critical
			{
				if(vState[ib] == 0){
					iFO = 1;
					vState[ib] = 1;
					rStartLoc = (long long)(ib) * (long long)(nReads2read);
					printf("%lld\n", rStartLoc);
					sprintf(outputFileLoc, "%s/data_%010lli.txt", outputFolder, rStartLoc + 1);
					bdLoc->fo = fopen(outputFileLoc, "w");

					bdLoc->indReadStart = rStartLoc;

					posL = 12 + ((long long)rStartLoc) * (long long)(readBlockSize);
					posR = 12 + ((long long)rStartLoc) * (long long)(readBlockSize);
#ifdef WIN32
					_fseeki64(fiLeft, posL, SEEK_SET);
					_fseeki64(fiRight, posR, SEEK_SET);
#else
					fseeko64(fiLeft, posL, SEEK_SET);
					fseeko64(fiRight, posR, SEEK_SET);
#endif
					fread(bdLoc->vmReadLeft, sizeof(__m128i), blockMsize * nReads2read, fiLeft);
					fread(bdLoc->vmReadRight, sizeof(__m128i), blockMsize * nReads2read, fiRight);
				}
			}

			if(iFO == 1) bdLoc->processReads(0, nReads2read);

#pragma omp critical
			{
				if(iFO == 1){
					fclose(bdLoc->fo); bdLoc->fo = nullptr;
				}
			}
		}
	}

	free(vState); vState = nullptr;

	if (nR > 0) {
		rStart = (long long)(nW) * (long long)(nReads2read);
		printf("%lli\n", rStart);
		sprintf(outputFile, "%s/data_%010lli.txt", outputFolder, rStart + 1);
		pBD[0]->fo = fopen(outputFile, "w");

		pBD[0]->indReadStart = rStart;

		posLL = 12 + ((long long)rStart) * (long long)(readBlockSize);
		posRR = 12 + ((long long)rStart) * (long long)(readBlockSize);

#ifdef WIN32
		_fseeki64(fiLeft, posLL, SEEK_SET);
		_fseeki64(fiRight, posRR, SEEK_SET);
#else
		fseeko64(fiLeft, posLL, SEEK_SET);
		fseeko64(fiRight, posRR, SEEK_SET);
#endif
		fread(pBD[0]->vmReadLeft, sizeof(__m128i), blockMsize * nR, fiLeft);
		fread(pBD[0]->vmReadRight, sizeof(__m128i), blockMsize * nR, fiRight);

		pBD[0]->processReads(0, nR);
		for(int i = 0; i < nR; i++){
			voBest[i + rStart] = pBD[0]->voBestLoc[i];
			voCountCandidates[i + rStart] = pBD[0]->voCountCandidatesLoc[i];
			voCountSolutions[i + rStart] = pBD[0]->voCountSolutionsLoc[i];
		}
		fclose(pBD[0]->fo); pBD[0]->fo = nullptr;
	}

	sprintf(outputFile, "%s/statInfo.txt", outputFolder);
	foStat = fopen(outputFile, "w");

	if (foStat == nullptr) {
		printf("Error: cannot open file \"%s\".\n", outputFile);
		return -9;
	}
	for (int i = 0; i < nReads; i++) {
		fprintf(foStat, "%i\t%i\t%i\t%i\n", i + 1, voCountCandidates[i], voCountSolutions[i], voBest[i]);
	}

	fclose(foStat); foStat = nullptr;

	return 0;
}

int main(int argc, char** argv) {
	int ires;
	contP* cp;
	cp = new contP();
	ires = cp->startProcessing(argc, argv);
	delete cp; cp = nullptr;

	return ires;
}

