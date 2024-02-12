//E:\Temp2\perlotSeeds\binaryBlocksText E:\Temp2\perlotSeeds\binaryBlocks 4 24 0 2

#include <fstream>
#include <string.h>
#include <chrono>
#include <smmintrin.h>
#include <emmintrin.h>

#include <omp.h>

//#define WIN32 1

const int nchunk = 1000000;

class SeedCheck {
public:
	SeedCheck(int nLen_in, int nMismatch_in, int* iv, int* ov);
	~SeedCheck();

	int findOrder();
	void checkCandidate(int istart);
	bool checkCandidatePrevious(int istart);

	long long count_Total;

	__m128i* uSeed, uRes[128], mUnit[256], mOO, mON, mZero, mOne;
	int* iOrder, * iOrderInv;
	int* inVec, * outVec, * ip;
	int nLen, nOne, nMismatch, nLQ;
	int nFound;
	bool isCheckPrevious;

	unsigned int u[4];

};

inline long long areSame(__m128i m1, __m128i m2) {
	/*
	__m128i m3, m4;
	m3 = _mm_xor_si128(m1, m2);
	m4 = _mm_bsrli_si128(m3, 8);
	m3 = _mm_or_si128(m3, m4);
	return _mm_extract_epi64(m3, 0);
	
	m3 = _mm_xor_si128(m1, m2);
	m4 = _mm_bsrli_si128(m3, 8);
	m3 = _mm_or_si128(m3, m4);*/

	
	//if nLen <= 64
	__m128i m3;
	m3 = _mm_xor_si128(m1, m2);

	return _mm_extract_epi64(m3, 0);
}

SeedCheck::SeedCheck(int nLen_in, int nMismatch_in, int* iv, int* ov) {
	nLen = nLen_in;
	nMismatch = nMismatch_in;

	count_Total = 0;

	nFound = 0;

	inVec = iv;
	outVec = ov;

	mOO = _mm_set1_epi32(0xffffffff);
	mZero = _mm_set1_epi32(0);

	for (int i = 0; i < nLen; i++) {
		u[0] = 0;
		u[1] = 0;
		u[2] = 0;
		u[3] = 0;

		u[i / 32] = 1 << (i % 32);
		mUnit[i] = _mm_set_epi32(u[3], u[2], u[1], u[0]);
		mUnit[i + nLen] = mUnit[i];
	}

	mOne = mZero;
	for (int i = 0; i < nLen; i++) {
		mOne = _mm_or_si128(mOne, mUnit[i]);
	}

	mON = mOO;
	for (int i = 0; i < nLen; i++) {
		mON = _mm_xor_si128(mON, mUnit[i]);
	}

	uSeed = (__m128i*)malloc(sizeof(__m128i) * nLen);
	ip = (int*)malloc(sizeof(int) * nLen);
	iOrder = (int*)malloc(sizeof(int) * nLen);
	iOrderInv = (int*)malloc(sizeof(int) * nLen);

	findOrder();

	nLQ = nLen - nMismatch + 1;
}

SeedCheck::~SeedCheck() {
	if (uSeed != nullptr) { free(uSeed); uSeed = nullptr; }
	if (ip != nullptr) { free(ip); ip = nullptr; }
	if (iOrder != nullptr) { free(iOrder); iOrder = nullptr; }
	if (iOrderInv != nullptr) { free(iOrderInv); iOrderInv = nullptr; }
}

int SeedCheck::findOrder() {
	int* vecT, * nCount, ibest, vbest, vq;

	vecT = (int*)malloc(sizeof(int) * 2 * nLen);
	nCount = (int*)malloc(sizeof(int) * nLen);
	iOrder[0] = 0;
	iOrderInv[0] = 0;

	for (int i = 0; i < 2 * nLen; i++) {
		vecT[i] = 0;
	}

	vecT[0] = 1;
	vecT[nLen] = 1;

	for (int i = 1; i < nLen; i++) {
		vbest = -1;
		ibest = -1;
		for (int j = 0; j < nLen; j++) {
			if (vecT[j] == 1) continue;
			vq = -1;
			for (int k = 1; k <= nLen; k++) {
				if (vecT[j + k] == 1) {
					vq = k;
					break;
				}
			}
			for (int k = 1; k <= nLen; k++) {
				if (vecT[nLen + j - k] == 1) {
					if (k < vq) vq = k;
					break;
				}
			}
			if (vq > vbest) {
				vbest = vq;
				ibest = j;
			}
		}
		vecT[ibest] = 1;
		vecT[ibest + nLen] = 1;
		iOrder[i] = ibest;
		iOrderInv[ibest] = i;
	}

	free(vecT); vecT = nullptr;
	free(nCount); nCount = nullptr;
	return 0;
}

void SeedCheck::checkCandidate(int istart) {
	int* ipos, ipnq, nq;
	bool isOk;

	__m128i m1, m2, m3, m4;

	ipos = inVec + nOne * istart;

	count_Total++;

	isCheckPrevious = true;

	for (int M = nMismatch; M <= nMismatch; M++) {
		nLQ = nLen - M + 1;
		m1 = mZero;
		for (int i = 0; i < nOne; i++) {
			m1 = _mm_or_si128(m1, mUnit[ipos[i]]);
		}
		uSeed[iOrderInv[0]] = m1;
		for (int i = 1; i < nLen; i++) {
			m2 = _mm_slli_epi64(m1, nLen - 1);
			m3 = _mm_srli_epi64(m1, 1);
			m4 = _mm_or_si128(m2, m3);
			m1 = _mm_and_si128(mOne, m4);
			uSeed[iOrderInv[i]] = m1;
		}

		isOk = true;

		for (int jj = 0; jj < M; jj++) {
			ip[jj] = jj;
		}

		uRes[0] = mZero;

		for (nq = 0; nq < M; ++nq) {
			uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ip[nq]]);
			if (areSame(uRes[nq + 1], mOne) == 0) {
				isOk = false;
				break;
			}
			if (areSame(uRes[nq + 1], uRes[nq]) == 0) break;
			if (nq == M - 1)break;
		}

		if (isOk) {
			ipnq = ip[nq];
			do {
				do {
					ipnq++;
					if (ipnq - nq == nLQ) {
						nq--;
						if (nq == 0) break;
						ipnq = ip[nq];
					}
					else {
						break;
					}
				} while (true);
				if (nq == 0) break;

				uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ipnq]);
				ip[nq] = ipnq;
				if (areSame(uRes[nq + 1], uRes[nq]) == 0) continue;

				while (nq < M - 1) {
					nq++;
					ipnq++;

					ip[nq] = ipnq;
					uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ipnq]);

					if (areSame(uRes[nq + 1], uRes[nq]) == 0)break;
				}

				if (areSame(uRes[nq + 1], mOne) == 0) {
					isOk = false;
					break;
				}
				if (!isOk) break;

			} while (true);
		}
		if (!isOk)break;

	}

	if (isOk) {
		memcpy(outVec + nFound * nOne, ipos, nOne * sizeof(int));
		nFound++;
	}
}

bool SeedCheck::checkCandidatePrevious(int istart) {
	int* ipos, ik, nq;

	ipos = inVec + nOne * istart;

	isCheckPrevious = true;

	nLQ = nLen - nMismatch + 1;

	for (int k = 0; k < nLen; k++) {
		uSeed[k] = mON;
		for (int i = 0; i < nOne; i++) {
			//ik = (ipos[i] + k) % nLen;
			ik = ipos[i] + k;
			uSeed[k] = _mm_or_si128(uSeed[k], mUnit[ik]);
		}
	}

	uRes[0] = mON;

	for (nq = 0; nq < nMismatch; ++nq) {

		uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ip[nq]]);
		if (areSame(uRes[nq + 1], mOO) == 0) return false;

		if (areSame(uRes[nq + 1], uRes[nq]) == 0) break;
		if (nq == nMismatch - 1)break;
	}

	return true;
}

void printP(FILE* fo, int* ipos, int nOne) {
	unsigned char qs[100];
	int ia, L, L8, iString[1000];

	L = ipos[nOne - 1] + 1;
	L8 = L / 8;
	if (L % 8 > 0) L8++;

	for (int i = 0; i < 8 * L8; i++) {
		iString[i] = 0;
	}

	for (int i = 0; i < nOne; i++) {
		iString[ipos[i]] = 1;
	}

	for (int i = 0; i < L8; i++) {
		ia = 0;

		for (int j = 7; j >= 0; j--) {
			ia = ia * 2 + iString[8 * i + j];
		}
		qs[i] = (unsigned char)ia;
	}
	fwrite(qs, sizeof(char), L8, fo);

	/*int m;
	m = 0;

	for (int i = 0; i < ipos[0]; i++) {
		fprintf(fo, "0");
	}
	fprintf(fo, "1");

	for (int i = 1; i < nOne; i++) {
		m = ipos[i] - ipos[i - 1] - 1;
		for (int j = 0; j < m; j++) {
			fprintf(fo, "0");
		}
		fprintf(fo, "1");
	}
	fprintf(fo, "\n");*/
}

int main2(char* inputFolder, char* outputFolder, int nMismatch, int nLen, int levelMin, int levelMax) {
	char fName[1000], sString[1000];

	int* nSolFound, ** pindSol, ** pcandSol;
	int* ipos, * gapSize, * indSame;
	int nL, nthreads, maxGap, mgg, nCount, nW;
	int nCut1, mgg2, nSame, nn, ndiff, ki, minWeight;

	long long countTest;

	bool tEqual, tLess, t, isFileClosed;

	SeedCheck** sc;

	FILE* fo, * fi;

	sprintf(fName, "%s/block_m_%i_len_%i_level_0.txt", inputFolder, nMismatch, nLen);
	fi = fopen(fName, "r");
	if (fi == nullptr) {
		printf("Error: cannot open  file \"%s\"\n", fName);
		return -1;
	}
	fscanf(fi, "%s", sString);
	fclose(fi); fi = nullptr;

	minWeight = 0;
	for (int i = 0; i < strlen(sString); i++) {
		if (sString[i] == '1') minWeight++;
	}

	printf("Weight (min): %i\n", minWeight);

#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}

	printf("Number of threads: %i\n", nthreads);

	sc = (SeedCheck**)malloc(sizeof(SeedCheck*) * nthreads);

	nSolFound = (int*)malloc(sizeof(int) * nthreads);
	pindSol = (int**)malloc(sizeof(int*) * nthreads);
	pcandSol = (int**)malloc(sizeof(int*) * nthreads);

	for (int i = 0; i < nthreads; i++) {
		pindSol[i] = (int*)malloc(sizeof(int) * nLen * nchunk);
		pcandSol[i] = (int*)malloc(sizeof(int) * nLen * nchunk);
	}

	for (int i = 0; i < nthreads; i++) {
		sc[i] = new SeedCheck(nLen, nMismatch, pcandSol[i], pindSol[i]);
	}


	ipos = (int*)malloc(sizeof(int) * nLen);
	gapSize = (int*)malloc(sizeof(int) * 2 * nLen);
	indSame = (int*)malloc(sizeof(int) * nLen);

	isFileClosed = true;
	fo = nullptr;

	for (int iLevel = levelMin; iLevel <= levelMax; iLevel++) {
		nW = minWeight - iLevel;

		printf("Mismatches: %i, len: %i, level: %i, weight: %i\n", nMismatch, nLen, iLevel, nW);
		if (nW < 1) continue;

		if (nW == 1) {
			nL = 0;
			ipos[0] = nLen - 1;

			sprintf(fName, "%s/block_m_%i_len_%i_level_%i.bin", outputFolder, nMismatch, nLen, iLevel);
			fo = fopen(fName, "wb");
			printP(fo, ipos, nW);
			fclose(fo);
			continue;
		}

		countTest = 0;
		nCount = 0;

		tEqual = false;
		tLess = false;

		nL = 0;
		ipos[nL] = nMismatch - 1;
		ipos[nW - 1] = nLen - 1;
		nCut1 = nLen - nW + 1;
		nn = (nW - 1) / 2;

		do {
			nCount = 0;
			do {
				do {
					ipos[nL]++;
					if (ipos[nL] - nL == nCut1) {
						nL--;
						if (nL < 0) break;
					}
					else {
						break;
					}
				} while (true);

				if (nL < 0) break;
				if (nL == 0) {
					maxGap = ipos[0] + 1;
					tEqual = false;
				}
				else {
					mgg = ipos[nL] - ipos[nL - 1];
					if (mgg > maxGap) {
						nL--;
						tEqual = false;
						tLess = true;
						continue;
					}
					if (mgg == maxGap) {
						tEqual = true;
					}
				}

				while (nL < nW - 2) {
					nL++;
					ipos[nL] = ipos[nL - 1] + 1;
				}

				t = false;
				for (int k = 0; k < nn; k++) {
					ndiff = (ipos[nW - 1 - k] - ipos[nW - 2 - k]) - (ipos[k + 1] - ipos[k]);
					if (ndiff == 0) continue;
					if (ndiff > 0) t = true;
					break;
				}

				if (t) continue;

				mgg2 = ipos[nW - 1] - ipos[nW - 2];

				if (mgg2 > maxGap) {
					tEqual = false;
					continue;
				}

				if (mgg2 == maxGap) {
					tEqual = true;
				}

				if (tEqual || tLess) {
					gapSize[0] = ipos[0] + 1;
					gapSize[nW] = gapSize[0];
					maxGap = gapSize[0];
					nSame = 1;
					indSame[0] = 0;
					t = false;

					for (int k = 1; k < nW; k++) {
						gapSize[k] = ipos[k] - ipos[k - 1];
						gapSize[k + nW] = gapSize[k];
						if (gapSize[k] == maxGap) {
							indSame[nSame] = k;
							nSame++;
						}
					}
					if (nSame == 1) {
						tLess = false;
						tEqual = false;
					}
					else {
						for (int k = 1; k < nSame; k++) {
							ki = indSame[k];
							for (int i = 1; i < nW; i++) {
								if (gapSize[ki + i] == gapSize[i])continue;
								if (gapSize[ki + i] > gapSize[i]) t = true;
								break;
							}
							if (t) break;
						}

						if (!t) {
							for (int k = 0; k < nSame; k++) {
								ki = indSame[k];
								for (int i = 0; i < nW; i++) {
									if (gapSize[ki - i + nW] == gapSize[i]) continue;
									if (gapSize[ki - i + nW] > gapSize[i]) t = true;
									break;
								}
								if (t) break;
							}
						}
					}

					if (t) continue;
				}

				for (int i = 0; i < nW; i++) {
					pcandSol[0][nCount * nW + i] = ipos[i];
				}

				nCount++;
				countTest++;

			} while (nCount < nchunk);

			printf("%lld\t%i\n", countTest, nCount);

			for (int i = 1; i < nthreads; i++) {
				memcpy(pcandSol[i], pcandSol[0], sizeof(int) * nCount * nW);
			}
			for (int i = 0; i < nthreads; i++) {
				sc[i]->nOne = nW;
				sc[i]->nFound = 0;
			}

#pragma omp parallel
			{
				int tid, i;
				SeedCheck* sc_loc;
				tid = omp_get_thread_num();
				sc_loc = sc[tid];
				sc_loc->isCheckPrevious = false;

#pragma omp for
				for (i = 0; i < nCount; i++) {
					sc_loc->checkCandidate(i);
				}
			}

			for (int m = 0; m < nthreads; m++) {
				if (sc[m]->nFound == 0) continue;

				for (int k = 0; k < sc[m]->nFound; k++) {
					if (isFileClosed) {
						sprintf(fName, "%s/block_m_%i_len_%i_level_%i.bin", outputFolder, nMismatch, nLen, iLevel);
						fo = fopen(fName, "wb");
						isFileClosed = false;
					}
					printP(fo, sc[m]->outVec + k * nW, nW);
				}
			}
		} while (nL >= 0);

		if (fo != nullptr) {
			fclose(fo); fo = nullptr;
			isFileClosed = true;
		}

		printf("Num: %lld\n", countTest);

		if (fo != nullptr) {
			fclose(fo); fo = nullptr;
		}

	}

	free(ipos); ipos = nullptr;
	free(nSolFound); nSolFound = nullptr;

	for (int i = 0; i < nthreads; i++) {
		free(pindSol[i]);
		free(pcandSol[i]);
		delete sc[i];
	}

	free(pindSol); pindSol = nullptr;
	free(pcandSol); pcandSol = nullptr;
	free(sc); sc = nullptr;

	return 0;
}

int printUsage() {
	printf("Arguments:\n");
	printf("1) Input folder\n");
	printf("2) Output folder\n");
	printf("3) Number of mismatches\n");
	printf("4) Size of blocks\n");
	printf("5) Minimum level\n");
	printf("6) Maximum level\n\n");
	return 0;
}

int main(int argc, char* argv[]) {
	int nMismatch, nLen, levelMin, levelMax;
	char outputFolder[1000], inputFolder[1000];
	int ires;

	if (argc != 7) {
		printUsage();
		return -99;
	}

	sprintf(inputFolder, "%s", argv[1]);
	sprintf(outputFolder, "%s", argv[2]);
	nMismatch = atoi(argv[3]);
	nLen = atoi(argv[4]);
	levelMin = atoi(argv[5]);
	levelMax = atoi(argv[6]);

	printf("1) Input folder: %s\n", inputFolder);
	printf("2) Output folder: %s\n", outputFolder);
	printf("3) Number of mismatches: %i\n", nMismatch);
	printf("4) Size of blocks: %i\n", nLen);
	printf("5) Minimum level: %i\n", levelMin);
	printf("6) Maximum level: %i\n\n", levelMax);

	auto start = std::chrono::steady_clock::now();

	ires = main2(inputFolder, outputFolder, nMismatch, nLen, levelMin, levelMax);

	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
#ifdef WIN32
	printf("It took me %lld ms.\n\n", elapsed.count());
#else
	printf("It took me %ld ms.\n\n", elapsed.count());
#endif

	return ires;
}
