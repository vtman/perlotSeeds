#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <string.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include <chrono>

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
//#define CASE_C68 1
//#define CASE_C72 1
//#define CASE_C76 1
//#define CASE_C80 1
//#define CASE_C84 1
//#define CASE_C88 1
//#define CASE_C92 1
//#define CASE_C96 1
//#define CASE_C100 1


//D:\Genome\Cram\ref38.acgt D:\Genome\library

const int nLetters = 4; //Modify
const int nBitsOut = (4 * nLetters);
const int nBytesOut = (nBitsOut / 8);
const int nFiles = (1 << (4 * nLetters));
const int indexLevel = 14;//Modify

const int inputBlockCount = 1000000;
const int outBlockPrint = 10000;


int printM(__m128i* vm, int n) {
	int* vi, q[4];
	vi = (int*)vm;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 32; j++) {
			for (int k = 0; k < 4; k++) {
				q[k] = (vi[4 * i + k] >> j) & 1;
			}
			if (q[0] == 1) {
				printf("A");
			}
			else if (q[1] == 1) {
				printf("C");
			}
			else if (q[2] == 1) {
				printf("G");
			}
			else if (q[3] == 1) {
				printf("T");
			}
			else {
				printf("_");
			}
		}
	}
	printf("\n");
	return 0;
}

int printOutput(__m128i* vm, int n) {
	int nWhole, nRem, * iv;
	int iAC, iAG, ij;

	iv = (int*)vm;
	nWhole = n / 32;
	nRem = n % 32;
	for (int i = 0; i < nWhole; i++) {
		for (int j = 0; j < 32; j++) {
			iAC = (iv[4 * i] >> j) & 1;
			iAG = (iv[4 * i + 1] >> j) & 1;
			if (iAC == 0 && iAG == 0)printf("T");
			if (iAC == 0 && iAG == 1)printf("G");
			if (iAC == 1 && iAG == 0)printf("C");
			if (iAC == 1 && iAG == 1)printf("A");
		}
	}

	for (int j = 0; j < nRem; j++) {
		iAC = (iv[4 * nWhole] >> j) & 1;
		ij = nRem + j;
		if (ij >= 32) {
			iAG = (iv[4 * nWhole + 1] >> (ij - 32)) & 1;
		}
		else {
			iAG = (iv[4 * nWhole] >> (ij)) & 1;
		}
		if (iAC == 0 && iAG == 0)printf("T");
		if (iAC == 0 && iAG == 1)printf("G");
		if (iAC == 1 && iAG == 0)printf("C");
		if (iAC == 1 && iAG == 1)printf("A");
	}

	printf("\n");
	return 0;
}

int printOriginal(__m128i* vm, int n) {
	int nWhole, nRem, * iv;
	int iA, iC, iG, iT;

	iv = (int*)vm;
	nWhole = n / 32;
	nRem = n % 32;
	for (int i = 0; i < nWhole; i++) {
		for (int j = 0; j < 32; j++) {
			iA = (iv[4 * i] >> j) & 1;
			iC = (iv[4 * i + 1] >> j) & 1;
			iG = (iv[4 * i + 2] >> j) & 1;
			iT = (iv[4 * i + 3] >> j) & 1;

			if (iA == 1)printf("A");
			if (iC == 1)printf("C");
			if (iG == 1)printf("G");
			if (iT == 1)printf("T");
			if (iA + iC + iG + iT == 0)printf("N");
		}
	}

	for (int j = 0; j < nRem; j++) {
		iA = (iv[4 * nWhole] >> j) & 1;
		iC = (iv[4 * nWhole + 1] >> j) & 1;
		iG = (iv[4 * nWhole + 2] >> j) & 1;
		iT = (iv[4 * nWhole + 3] >> j) & 1;

		if (iA == 1)printf("A");
		if (iC == 1)printf("C");
		if (iG == 1)printf("G");
		if (iT == 1)printf("T");
		if (iA + iC + iG + iT == 0)printf("N");
	}

	printf("\n");
	return 0;
}


int printData(char* vD, int n) {
	for (int i = 0; i < n; i++) {
		printf("%1x%1x", vD[i] & 15, (vD[i] >> 4) & 15);
	}
	return 0;
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

#ifdef CASE_C68
	//C68
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x0000000f);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 4));
#endif

#ifdef CASE_C72
	//C72
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x000000ff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 8));
#endif

#ifdef CASE_C76
	//C76
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x00000fff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 12));
#endif

#ifdef CASE_C80
	//C80
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x0000ffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 16));
#endif

#ifdef CASE_C84
	//C84
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x000fffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 20));
#endif

#ifdef CASE_C88
	//C88
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x00ffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 24));
#endif

#ifdef CASE_C92
	//C92
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mMask = _mm_set_epi32(0, 0, 0, 0x0fffffff);
	mAC = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4)));
	mAG = _mm_and_si128(mMask, _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8)));
	vPattern[1] = _mm_or_si128(mAC, _mm_slli_epi64(mAG, 28));
#endif

#ifdef CASE_C96
	//C96
	__m128i mAC, mAG, mMask, res[2];

	mAC = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 4));
	mAG = _mm_or_si128(m[0], _mm_bsrli_si128(m[0], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);

	mAC = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 4));
	mAG = _mm_or_si128(m[1], _mm_bsrli_si128(m[1], 8));
	res[1] = _mm_unpacklo_epi32(mAC, mAG);

	vPattern[0] = _mm_unpacklo_epi64(res[0], res[1]);

	mAC = _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 4));
	mAG = _mm_or_si128(m[2], _mm_bsrli_si128(m[2], 8));
	res[0] = _mm_unpacklo_epi32(mAC, mAG);
	vPattern[1] = _mm_unpacklo_epi64(res[0], _mm_set1_epi32(0));
#endif



	return 0;
}


int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], outputFolder[1000], outputTemplate[1000], mString[1000];
	char sSeed[1000];
	char* vQ;
	char** pcData;

	int nSeed, nBitsPerRecord, nBytesPerRecord;
	int blockStep, nWhole, nRem;
	int n8, n16, n32, n128;
	int uBlockMask, uFile, uStart;
	int recordsPerBlock, memoryPerBlock, nBlocks;
	int nWeight1, nWeight2, nDoubleWeight;
	int nBP, AA;

	int mOne[32];

	int* iFileState, * ivCountRecords;

	unsigned int uPos;

	FILE* fo, * fi;

	__m128i mRef[10], mOut[10], vPattern[10];
	__m128i* vi;
	__m128i mZero, m1, m2, m3, m4, m5;

	bool t;

	long long fDiff, fPos, uC, fileSize;

	auto start = std::chrono::steady_clock::now();

	if (argc != 3) {
		printf("Error: wrong number of input parameters.\n");
		printf("1) path to ACGT file\n");
		printf("2) output folder\n");
		return -1;
	}

	iFileState = (int*)malloc(sizeof(int) * nFiles);

#ifdef CASE_T0V1
	sprintf(sSeed, "######_#######_#######_#######_#######_#######_#######_#######_######");//T0V1
#endif
#ifdef CASE_T0V2
	sprintf(sSeed, "#_#####_###__#_#####_###__#_#####_###__#_#####_###__#_#####_###");//T0V2
#endif
#ifdef CASE_T0V3
	sprintf(sSeed, "##_#_####___#__##_#_####___#__##_#_####___#__##_#_####___#__##");//T0V3
#endif

#ifdef CASE_T1V0
	sprintf(sSeed, "######@#######@#######@#######@#######@#######@#######@#######@######");//T1V0
#endif
#ifdef CASE_T1V1
	sprintf(sSeed, "#####_@##@#@#####_@##@#@#####_@##@#@#####_@##@#@#####_@##@#@#####");//T1V1
#endif
#ifdef CASE_T1V2
	sprintf(sSeed, "##_#_####__@#@_##_#_####__@#@_##_#_####__@#@_##_#_####__@#@_##");//T1V2
#endif
#ifdef CASE_T1V3
	sprintf(sSeed, "##_@_@#@#___#__##_@_@#@#___#__##_@_@#@#___#__##_@_@#@#___#__##");//T1V3
#endif

#ifdef CASE_T2V0
	sprintf(sSeed, "#@#####@###@@#@#####@###@@#@#####@###@@#@#####@###@@#@#####@###");//T2V0
#endif
#ifdef CASE_T2V1
	sprintf(sSeed, "##@#@####_@@#@@##@#@####_@@#@@##@#@####_@@#@@##@#@####_@@#@@##");//T2V1
#endif
#ifdef CASE_T2V2
	sprintf(sSeed, "##@#@@@##_#_##@@__@@_##@#@@@##_#_##@@__@@_##@#@@@##_#_##");//T2V2
#endif

#ifdef CASE_T3V0
	sprintf(sSeed, "##@#@####@@@#@@##@#@####@@@#@@##@#@####@@@#@@##@#@####@@@#@@##");//T3V0
#endif
#ifdef CASE_T3V1
	sprintf(sSeed, "##@@@@#_@@##@@@@#_@@##@@@@#_@@##@@@@#_@@##@@@@#_@@##@@@@#_@@##@@@@#");//T3V1
#endif
#ifdef CASE_T3V2
	sprintf(sSeed, "#@@##@__@#@_#_@@_@#@@##@__@#@_#_@@_@#@@##@__@#@_#_@@_@#@@##");//T3V2
#endif

#ifdef CASE_T4V0
	sprintf(sSeed, "#@#####@##@@###@@@@##@#@#@@#@@@#@#####@##@@###");//T4V0
#endif
#ifdef CASE_T4V1
	sprintf(sSeed, "#@@@#@#@#_@@@@#@@@#@#@#_@@@@#@@@#@#@#_@@@@#@@@#@#@#_@@@@#@@@#@#");//T4V1
#endif

#ifdef CASE_T5V0
	sprintf(sSeed, "##@@@##@#@#####@@@@@@@@@@@@@@@@##@@@##@#@#####");//T5V0
#endif
#ifdef CASE_T5V1
	sprintf(sSeed, "#@@@@@##_@@@@@#@@@@@##_@@@@@#@@@@@##_@@@@@#@@@@@##_@@@@@#@@@@@#");//T5V1
#endif

#ifdef CASE_C16
	sprintf(sSeed, "################");//B16
#endif
#ifdef CASE_C20
	sprintf(sSeed, "####################");//B20
#endif
#ifdef CASE_C24
	sprintf(sSeed, "########################");//B24
#endif
#ifdef CASE_C28
	sprintf(sSeed, "############################");//B28
#endif
#ifdef CASE_C32
	sprintf(sSeed, "################################");//B32
#endif
#ifdef CASE_C36
	sprintf(sSeed, "####################################");//B36
#endif
#ifdef CASE_C40
	sprintf(sSeed, "########################################");//B40
#endif
#ifdef CASE_C44
	sprintf(sSeed, "############################################");//B44
#endif
#ifdef CASE_C48
	sprintf(sSeed, "################################################");//B48
#endif
#ifdef CASE_C52
	sprintf(sSeed, "####################################################");//B52
#endif
#ifdef CASE_C56
	sprintf(sSeed, "########################################################");//B56
#endif
#ifdef CASE_C60
	sprintf(sSeed, "############################################################");//B60
#endif
#ifdef CASE_C64
	sprintf(sSeed, "################################################################");//B64
#endif
#ifdef CASE_C68
	sprintf(sSeed, "####################################################################");//B68
#endif
#ifdef CASE_C72
	sprintf(sSeed, "########################################################################");//B72
#endif
#ifdef CASE_C76
	sprintf(sSeed, "############################################################################");//B76
#endif
#ifdef CASE_C80
	sprintf(sSeed, "################################################################################");//B80
#endif
#ifdef CASE_C84
	sprintf(sSeed, "####################################################################################");//B84
#endif
#ifdef CASE_C88
	sprintf(sSeed, "########################################################################################");//B88
#endif
#ifdef CASE_C92
	sprintf(sSeed, "############################################################################################");//B92
#endif
#ifdef CASE_C96
	sprintf(sSeed, "################################################################################################");//B96
#endif



	sprintf(inputFile, "%s", argv[1]);
	sprintf(outputFolder, "%s", argv[2]);
	nSeed = strlen(sSeed);

	printf("Input file: \"%s\".\n", inputFile);
	printf("Output folder: \"%s\".\n", outputFolder);
	printf("Length of a seed: %i\n\n", nSeed);

	nWeight1 = 0;
	nWeight2 = 0;
	for (int i = 0; i < nSeed; i++) {
		if (sSeed[i] == '#') nWeight1++;
		if (sSeed[i] == '@') nWeight2++;
	}

	nDoubleWeight = 2 * nWeight1 + nWeight2;

	sprintf(outputTemplate, "%s/original/%%0%ix.bin", outputFolder, nLetters);
	printf("Output template: %s\n", outputTemplate);

	sprintf(mString, "%s/info.txt", outputFolder);

	fo = fopen(mString, "w");
	if (fo == nullptr) {
		printf("Error: cannot open file %s.\n", mString);
		return -1;
	}
	fprintf(fo, "SeedLength\t%i\n", nSeed);
	fprintf(fo, "DoubleWeight\t%i\n", nDoubleWeight);
	fprintf(fo, "Seed\t");
	fprintf(fo, "%s", sSeed);
	fprintf(fo, "\n");
	fprintf(fo, "Letters\t%i\n", nLetters);
	fprintf(fo, "IndexLevel\t%i\n", indexLevel);
	fclose(fo); fo = nullptr;

	nBitsPerRecord = 32 + nDoubleWeight;

	nBP = nDoubleWeight / 8;
	if (nDoubleWeight % 8 > 0) nBP++;

	printf("Bits per record: %i\n", nBitsPerRecord);

	nBytesPerRecord = nBitsPerRecord - nBitsOut;
	if (nBytesPerRecord % 8 > 0) nBytesPerRecord += 8;
	nBytesPerRecord /= 8;

	printf("Bytes per record (storage): %i\n", nBytesPerRecord);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\".\n", inputFile);
		return -2;
	}

#ifdef WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
#else
	fseeko64(fi, 0, SEEK_END);
	fileSize = ftello64(fi);
#endif

	printf("File size: %lli\n", fileSize);

	recordsPerBlock = 32 * inputBlockCount;
	printf("Records per block: %i\n", recordsPerBlock);

	n16 = (nSeed - 1) / 32;
	if ((nSeed - 1) % 32 > 0) n16++;
	memoryPerBlock = 16 * n16;
	memoryPerBlock += 16 * inputBlockCount;
	printf("Memory per block: %i\n", memoryPerBlock);

	blockStep = 16 * inputBlockCount;
	printf("Block step: %i\n", blockStep);

	nBlocks = (int)(fileSize / ((long long)blockStep));
	if (((long long)nBlocks) * ((long long)blockStep) < fileSize) nBlocks++;

	printf("Number of blocks: %i\n", nBlocks);

	n32 = inputBlockCount;

	n128 = nSeed;
	if (n128 % 32 > 0) n128 += 32;
	n128 /= 32;

#ifdef WIN32
	vQ = (char*)_aligned_malloc(sizeof(__m128i) * n128, sizeof(__m128i));
#else
	vQ = (char*)aligned_alloc(sizeof(__m128i), sizeof(__m128i) * n128);
#endif

	n8 = nDoubleWeight;
	if (n8 % 8 > 0) n8 += 8;
	n8 /= 8;


	for (int i = 0; i < 32; i++) {
		mOne[i] = 1 << i;
	}

	uBlockMask = 0;
	for (int i = 0; i < nBitsOut; i++) {
		uBlockMask |= mOne[i];
	}

	nWhole = nSeed / 32;
	nRem = nSeed % 32;

	if (nRem > 0) nWhole++;

	for (int i = 0; i < nWhole; i++) {
		AA = 0;
		for (int j = 0; j < 32; j++) {
			if (32 * i + j >= nSeed) break;
			if (sSeed[32 * i + j] == '_')continue;
			AA |= (1 << j);
		}
		mRef[i] = _mm_set1_epi32(AA);
	}

	mZero = _mm_set1_epi32(0);

	uC = 0;

#ifdef WIN32
	vi = (__m128i*)_aligned_malloc((n32 + n16) * sizeof(__m128i), sizeof(__m128i));
#else
	vi = (__m128i*)aligned_alloc(sizeof(__m128i), (n32 + n16) * sizeof(__m128i));
#endif

	uPos = 0;

	iFileState = (int*)malloc(sizeof(int) * nFiles);
	ivCountRecords = (int*)malloc(sizeof(int) * nFiles);
	pcData = (char**)malloc(sizeof(char*) * nFiles);

	for (int i = 0; i < nFiles; i++) {
		iFileState[i] = -1;
		ivCountRecords[i] = 0;
		pcData[i] = (char*)malloc(sizeof(char) * nBytesPerRecord * outBlockPrint);
	}

	for (int ib = 0; ib < nBlocks; ib++) {
		fPos = ((long long)ib) * ((long long)blockStep);

#ifdef WIN32
		_fseeki64(fi, fPos, SEEK_SET);
#else
		fseeko64(fi, fPos, SEEK_SET);
#endif

		fDiff = fileSize - fPos;
		if (fDiff <= 0)break;
		if (fDiff < (long long)memoryPerBlock) {
			n32 = (int)(fDiff / 16) - n16;
			memoryPerBlock = n32 * 16;
		}
		fread(vi, sizeof(__m128i), n32 + n16, fi);
		//printf("Size m128: %i\n", sizeof(__m128i));

		for (int ii = 0; ii < n32; ii++) {
			for (int j = 0; j < 32; j++) {
				m3 = mZero;
				for (int k = 0; k < nWhole; k++) {
					m1 = _mm_srli_epi32(vi[ii + k], j);
					m2 = _mm_slli_epi32(vi[ii + k + 1], 32 - j);
					mOut[k] = _mm_and_si128(mRef[k], _mm_or_si128(m1, m2));

					m4 = _mm_or_si128(mOut[k], _mm_bsrli_si128(mOut[k], 8));
					m5 = _mm_or_si128(m4, _mm_bsrli_si128(m4, 4));
					m3 = _mm_or_si128(m3, _mm_xor_si128(m5, mRef[k]));
				}

				formSignatureM(mOut, vPattern);

				for (int k = 0; k < nBP; k++) {
					vQ[k] = ((char*)vPattern)[k];
				}

				t = (_mm_extract_epi32(m3, 0) == 0);

				if (t) {
					uFile = uBlockMask & (*(int*)vQ);
					//printf("%i\t%08x\t", j, uFile);
					uStart = ivCountRecords[uFile] * nBytesPerRecord;
					for (int i = 0; i < n8 - nBytesOut; i++) {
						pcData[uFile][uStart + i] = vQ[nBytesOut + i];
					}
					for (int i = 0; i < 4; i++) {
						pcData[uFile][uStart + n8 - nBytesOut + i] = *(((char*)&uPos) + i);
					}
					ivCountRecords[uFile]++;
					if (ivCountRecords[uFile] == outBlockPrint) {
						sprintf(outputFile, outputTemplate, uFile);
						if (iFileState[uFile] == -1) {
							fo = fopen(outputFile, "wb");
							iFileState[uFile] = 0;
						}
						else {
							fo = fopen(outputFile, "ab");
						}
						fwrite(pcData[uFile], sizeof(char), ivCountRecords[uFile] * nBytesPerRecord, fo);
						ivCountRecords[uFile] = 0;
						fclose(fo); fo = nullptr;
					}
				}
				uC++;
				uPos++;
			}

		}

		printf("%i\t%i\n", ib, n32);
	}

	for (int i = 0; i < nFiles; i++) {
		sprintf(outputFile, outputTemplate, i);
		if (iFileState[i] == -1) {
			fo = fopen(outputFile, "wb");
			iFileState[i] = 0;
		}
		else {
			fo = fopen(outputFile, "ab");
		}
		fwrite(pcData[i], sizeof(char), ivCountRecords[i] * nBytesPerRecord, fo);
		ivCountRecords[i] = 0;
		fclose(fo); fo = nullptr;
	}

	for (int i = 0; i < nFiles; i++) {
		free(pcData[i]); pcData[i] = nullptr;
	}
	free(pcData); pcData = nullptr;
	free(ivCountRecords); ivCountRecords = nullptr;
	free(iFileState); iFileState = nullptr;

	printf("Total u: %lli\n", uC);

#ifdef WIN32
	_aligned_free(vi); vi = nullptr;
	_aligned_free(vQ); vQ = nullptr;
#else
	free(vi); vi = nullptr;
	free(vQ); vQ = nullptr;
#endif

	fclose(fi); fi = nullptr;

	free(iFileState); iFileState = nullptr;

	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
#ifdef WIN32
	printf("It took me %lld ms.\n\n", elapsed.count());
#else
	printf("It took me %ld ms.\n\n", elapsed.count());
#endif

	return 0;
}
