#define _CRT_SECURE_NO_WARNINGS

#define DSIGMA 4
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define FALSE	0
#define TRUE	1
#define XSIZE	4200
#define WSIZE	256
#define SIGMA	256
#define UNDEFINED	-1
#define HALFDEFINED	-2
#define WORD	32
#define OUTPUT(j)	count++

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <chrono>

using namespace std;

class findExactWithMAS
{
public:
	findExactWithMAS()
	{
		containerString = "";
		readRefDNA.open("refDNA.txt");
		readShortReads.open("shortReadList.txt");
		readRefDNA >> refDNA;
		referenceDNASize = refDNA.length();
		for (int i = 0; i < referenceDNASize; i++)
		{
			containerString += 'N';
		}
	}
	~findExactWithMAS()
	{
		readRefDNA.close();
		readShortReads.close();
	}

	int MAS(string Text, string Pattern, int start);
	void divideAndCompare(string shortRead);
	bool checkOtherStrings(int index, int position, string shortRead);
	void run();
	void compareWithRef();
	void saveContainer();

private:
	string refDNA;
	string shortReads;
	string* shortReadsArray;
	string containerString;

	int shortReadSize;
	int dividedShortReadSize;
	int snpPercentage;
	int misMatchNum;
	int referenceDNASize;
	int i, l, k, s, w, nMinusm;
	int max_pos, avr_shift, max_avr_shift;
	bool match;
	int count;
	int scan[XSIZE], shift[XSIZE][SIGMA];
	int U[XSIZE], safe[XSIZE + 1];

	ifstream readRefDNA;
	ifstream readShortReads;
};

int findExactWithMAS::MAS(string T, string P, int start)
{
	// 전처리 파트

	unsigned char* Pattern = (unsigned char*)P.c_str();
	unsigned char* Text = (unsigned char*)T.c_str();

	char dna[DSIGMA] = { 'A', 'C', 'G', 'T' };

	int freq[SIGMA];
	freq['A'] = 293;	freq['C'] = 207;
	freq['G'] = 207;	freq['T'] = 293;

	// U 초기화
	for (l = 0; l < dividedShortReadSize; l++)
	{
		U[l] = 1;
		// shift 배열 초기화, DSIGMA == 4
		for (s = 0; s < DSIGMA; s++)
			shift[l][dna[s]] = 1; // 아스키 코드로 A C G T 배열 인덱스 사용
	}

	// safe 배열 초기화, 패턴 길이만큼, 
	for (k = 1; k <= dividedShortReadSize; k++)
		safe[k] = 0;

	for (i = 0; i < dividedShortReadSize; i++)
	{
		for (l = 0; l < dividedShortReadSize; l++)
		{
			if (U[l] == 1) {
				// A C G T 차례대로 
				for (s = 0; s < DSIGMA; s++)
				{
					// 패턴과 비교하여 shift[l][dna[s]] 값 계산하는 반복문
					for (k = shift[l][dna[s]]; k <= dividedShortReadSize; k++)
					{
						if (safe[k] == 0 && ((l - k < 0) || dna[s] == (Pattern[l - k])))
						{
							shift[l][dna[s]] = k;
							break;
						}
					}
				}
			}
		}

		// max average shift 계산 반복문 
		max_avr_shift = 0;
		for (l = 0; l < dividedShortReadSize; l++)
		{
			if (U[l] == 1) {
				avr_shift = 0;
				// average shift 먼저 계산 후
				for (s = 0; s < DSIGMA; s++) {
					avr_shift = avr_shift + shift[l][dna[s]] * freq[dna[s]];

				} // 값이 가장 큰 max average shift 값 구하기
				if ((max_avr_shift < avr_shift) || (max_avr_shift == avr_shift && freq[Pattern[max_pos]] > freq[Pattern[l]])) {
					max_avr_shift = avr_shift;
					max_pos = l;
				}
			}
		}

		scan[i] = max_pos;
		U[max_pos] = 0;

		for (k = 1; k <= max_pos; k++)
		{
			if (Pattern[max_pos] != Pattern[max_pos - k])
			{
				safe[k] = 1;
			}
		}
	}

	memcpy(Text + T.length(), Pattern, P.length());
	w = start;
	nMinusm = referenceDNASize - dividedShortReadSize;
	while (1)
	{
		while (Pattern[(l = scan[0])] != (s = Text[w + l])) {
			w += shift[l][s];
		}
		if (w <= nMinusm) {
			i = 1;
			while (i < dividedShortReadSize && Pattern[(l = scan[i])] == (s = Text[w + l])) {
				i++;
			}
			if (i == dividedShortReadSize) {
				memset(Pattern, 0, sizeof(Pattern));
				memset(Text, 0, sizeof(Text));
				return w;
			}
			w += shift[l][s];
		}
		else {
			memset(Pattern, 0, sizeof(Pattern));
			memset(Text, 0, sizeof(Text));
			return -1;
		}
	}
}

void findExactWithMAS::run()
{
	string inLine;
	int position, i = 0;
	cout << "숏리드의 길이와 스닙의 확률을 입력해주세요 :: ";
	cin >> shortReadSize >> snpPercentage;

	dividedShortReadSize = shortReadSize / (snpPercentage + 1);
	misMatchNum = (shortReadSize / 100) * (snpPercentage) * 2;

	cout << "분석 시작" << endl;

	chrono::system_clock::time_point StartTime01 = chrono::system_clock::now();
	while (getline(readShortReads, inLine))
	{
		cout << i << endl;
		divideAndCompare(inLine);
		i++;
	}

	chrono::system_clock::time_point EndTime01 = chrono::system_clock::now();
	chrono::milliseconds milli01 = std::chrono::duration_cast<chrono::milliseconds> (EndTime01 - StartTime01);
	cout << endl << "MAS Time Count L " << milli01.count() << "milli seconds" << endl;
	cout << "숏리드의 길이 : " << shortReadSize << endl;
	cout << "읽은 숏리드의 갯수 : " << i << "개" << endl;
	cout << "허용 mismatch 수 : " << misMatchNum << endl;
}

void findExactWithMAS::divideAndCompare(string shortRead) {
	int position;
	shortReadsArray = new string[snpPercentage + 1];
	match = false;

	for (int i = 0; i < snpPercentage + 1; i++)
	{
		shortReadsArray[i] = shortRead.substr(i*dividedShortReadSize, dividedShortReadSize);
	}

	for (int i = 0; i < snpPercentage + 1; i++)
	{
		int patternSize = dividedShortReadSize;
		int textSize = referenceDNASize;

		int j = 0;

		int index = this->MAS(refDNA, shortReadsArray[i], 0);

		if (index != -1)
		{
			while (!checkOtherStrings(index, i, shortRead))
			{
				index = this->MAS(refDNA, shortReadsArray[i], index + 1);
				if (index == -1) break;
			}
			if (match)
			{
				break;
			}
			else
				continue;
		}
	}

	delete[] shortReadsArray;
}

bool findExactWithMAS::checkOtherStrings(int index, int position, string shortRead)
{
	int compare = (index - (position * dividedShortReadSize));
	int misMatchCount = 0;
	if (index == -1)
		return false;
	else if (compare > 0)
	{
		for (int i = 0; i < shortRead.length(); i++)
		{
			if (refDNA[compare + i] != shortRead[i])
				misMatchCount += 1;
			if (misMatchCount > misMatchNum)
				break;
		}
		if (misMatchCount > misMatchNum)
			return false;
		else {
			for (int i = 0; i < shortRead.length(); i++)
			{
				containerString[compare + i] = shortRead[i];
			}
			match = true;
			return true;
		}
	}
}

void findExactWithMAS::saveContainer()
{
	ofstream writeConstructed("constructed.txt");
	writeConstructed << containerString;
	writeConstructed.close();
}

void findExactWithMAS::compareWithRef()
{
	ifstream mySequenceDNA("myDNA.txt");
	string comp;
	mySequenceDNA >> comp;

	int cnt01 = 0, cnt02 = 0;
	for (int i = 0; i < refDNA.length(); i++)
	{
		if (refDNA[i] != containerString[i])
			cnt01 += 1;
		if (comp[i] != containerString[i])
			cnt02 += 1;
	}

	cout << "reference DNA length : " << referenceDNASize << endl;
	cout << "reference mismatch : " << cnt01 << endl;
	cout << "my Sequence mismatch : " << cnt02 << endl;
	mySequenceDNA.close();
}

int main()
{
	findExactWithMAS* aa = new findExactWithMAS();
	clock_t clockStart = clock();
	aa->run();
	aa->compareWithRef();
	double dueTime = ((float)clock() - clockStart) / CLOCKS_PER_SEC;
	double sec = dueTime;
	cout << "걸린 시간 :: " << sec << "초\n";
	aa->saveContainer();
}