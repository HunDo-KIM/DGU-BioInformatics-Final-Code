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
	// ��ó�� ��Ʈ

	unsigned char* Pattern = (unsigned char*)P.c_str();
	unsigned char* Text = (unsigned char*)T.c_str();

	char dna[DSIGMA] = { 'A', 'C', 'G', 'T' };

	int freq[SIGMA];
	freq['A'] = 293;	freq['C'] = 207;
	freq['G'] = 207;	freq['T'] = 293;

	// U �ʱ�ȭ
	for (l = 0; l < dividedShortReadSize; l++)
	{
		U[l] = 1;
		// shift �迭 �ʱ�ȭ, DSIGMA == 4
		for (s = 0; s < DSIGMA; s++)
			shift[l][dna[s]] = 1; // �ƽ�Ű �ڵ�� A C G T �迭 �ε��� ���
	}

	// safe �迭 �ʱ�ȭ, ���� ���̸�ŭ, 
	for (k = 1; k <= dividedShortReadSize; k++)
		safe[k] = 0;

	for (i = 0; i < dividedShortReadSize; i++)
	{
		for (l = 0; l < dividedShortReadSize; l++)
		{
			if (U[l] == 1) {
				// A C G T ���ʴ�� 
				for (s = 0; s < DSIGMA; s++)
				{
					// ���ϰ� ���Ͽ� shift[l][dna[s]] �� ����ϴ� �ݺ���
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

		// max average shift ��� �ݺ��� 
		max_avr_shift = 0;
		for (l = 0; l < dividedShortReadSize; l++)
		{
			if (U[l] == 1) {
				avr_shift = 0;
				// average shift ���� ��� ��
				for (s = 0; s < DSIGMA; s++) {
					avr_shift = avr_shift + shift[l][dna[s]] * freq[dna[s]];

				} // ���� ���� ū max average shift �� ���ϱ�
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
	cout << "�������� ���̿� ������ Ȯ���� �Է����ּ��� :: ";
	cin >> shortReadSize >> snpPercentage;

	dividedShortReadSize = shortReadSize / (snpPercentage + 1);
	misMatchNum = (shortReadSize / 100) * (snpPercentage) * 2;

	cout << "�м� ����" << endl;

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
	cout << "�������� ���� : " << shortReadSize << endl;
	cout << "���� �������� ���� : " << i << "��" << endl;
	cout << "��� mismatch �� : " << misMatchNum << endl;
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
	cout << "�ɸ� �ð� :: " << sec << "��\n";
	aa->saveContainer();
}