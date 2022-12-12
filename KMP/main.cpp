#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <chrono>

using namespace std;

class findExactWithKMP
{
public:
	findExactWithKMP()
	{
		containerString = "";
		readRefDNA.open("refDNA.txt");
		readShortReads.open("shortReadList.txt");
		readRefDNA >> refDNA;
		referenceDNASize = refDNA.length();
		for (int i = 0; i < refDNA.length(); i++)
		{
			containerString += 'N';
		}
	}
	~findExactWithKMP()
	{
		readRefDNA.close();
		readShortReads.close();
	}
	int KMP(string T, string P, int start, vector<int> SP);

	void divideAndCompare(string shortRead);
	bool checkOtherStrings(int index, int position, string shortRead);
	void run();
	void compareWithReference();
	void saveContainer();

private:
	string refDNA;	// 레퍼런스 DNA
	string shortReads; // 숏리드
	string* shortReadsArray; // 나눠진 숏리드를 담는 배열
	string containerString;
	vector<int> SP;

	int shortReadSize;
	int dividedShortReadSize;
	int snpPercentage;
	int misMatchNum;
	int referenceDNASize;

	bool match;

	ifstream readRefDNA;
	ifstream readShortReads;
};


void findExactWithKMP::run()
{
	string inLine;
	int position, i = 0;

	cout << "숏리드의 길이와 스닙의 확률(%)을 입력해주세요 :: ";
	cin >> shortReadSize >> snpPercentage;
	dividedShortReadSize = shortReadSize / (snpPercentage + 1); // 숏리드 쪼개기
	misMatchNum = (shortReadSize / 100) * (snpPercentage) * 2; // 허용 mismatch 수
	cout << "분석 시작 " << endl;
	chrono::system_clock::time_point StartTime01 = chrono::system_clock::now();
	while (getline(readShortReads, inLine))
	{
		cout << i << endl;
		divideAndCompare(inLine); // KMP 돌리는 중
		i++;
	}
	chrono::system_clock::time_point EndTime01 = chrono::system_clock::now();
	chrono::milliseconds milli01 = std::chrono::duration_cast<chrono::milliseconds>(EndTime01 - StartTime01);
	cout << endl << "KMP Time Count : " << milli01.count() << "milli second" << endl;
	cout << "숏리드의 길이 : " << shortReadSize;
	cout << "읽은 숏리드의 갯수 : " << i << "개" << endl;
	cout << "허용 mismatch 수 : " << misMatchNum << endl;
}


void findExactWithKMP::divideAndCompare(string shortRead) {
	int position; // 나눠진 숏리드가 원본 숏리드의 어떤 위치에 있는지

	shortReadsArray = new string[(snpPercentage + 1)]; // 나눠진 숏리드 문자열을 담을 배열.
	match = false;

	// 원본 숏리드를 나눠 담는 반복문
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		shortReadsArray[i] = shortRead.substr(i * dividedShortReadSize, dividedShortReadSize);
	}

	// 나눠진 숏리드들을 하나하나 검색함. 퍼펙트 매칭을 구하는 과정임.
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		SP.clear();
		// 텍스트의 길이와 패턴의 길이
		int patternSize = dividedShortReadSize;
		int textSize = referenceDNASize;
		// 패턴의 길이 만큼 배열을 생성
		SP.assign(patternSize, 0);
		int j = 0;

		// i는 1부터 시작, j는 0부터 시작
		for (int k = 1; k < patternSize; k++)
		{
			// 만약 j가 0보다 크면서, 패턴의 i위치와 j위치가 다를 경우 j를 뒤로 하나씩 이동
			while (j > 0 && shortReadsArray[i][k] != shortReadsArray[i][j])
			{
				j = SP[j - 1];
			}

			// j는 1증가, i위치에 추가.
			if (shortReadsArray[i][k] == shortReadsArray[i][j])
			{
				SP[k] = ++j;
			}
		}

		// SP 생성 완료
		j = 0;

		int index = this->KMP(refDNA, shortReadsArray[i], 0, SP);

		if (index != -1)
		{
			while (!checkOtherStrings(index, i, shortRead))
			{
				index = this->KMP(refDNA, shortReadsArray[i], index + 1, SP);
				// 다른 인덱스 검색
				if (index == -1) break; // 검색 못함.
			}
			if (match)
			{
				// 인덱스를 발견했을 때
				break;
			}
			else
				continue;
		}
	}

	delete[] shortReadsArray;

}


int findExactWithKMP::KMP(string Text, string Pattern, int start, vector<int> SP)
{
	int j = 0;
	// 0 부터 텍스트 길이 끝까지
	for (int i = start; i < referenceDNASize; i++)
	{
		// 접두사 접미사 테이블의 j 인덱스의 글자가 텍스트의 i 인덱스의 글자와 일치할 때 까지.
		while (j > 0 && Text[i] != Pattern[j])
		{
			j = SP[j - 1];
		}

		// 만약 글자가 같을 경우
		if (Text[i] == Pattern[j])
		{
			// 글자가 모두 일치할 경우
			if (j == dividedShortReadSize - 1)
			{
				return i - dividedShortReadSize + 1;
			}
			// 글자가 일치하는 경우 패턴 끝까지 비교함, j++
			else
			{
				j++;
			}
		}
	}
	return -1;
}

bool findExactWithKMP::checkOtherStrings(int index, int position, string shortRead)
{
	int compare = (index - (position * dividedShortReadSize));
	int misMatchCount = 0;
	if (index == -1)
		return false;
	else if (compare > 0) // 해당 위치에서 숏리드를 비교할 수 있는지
	{
		for (int i = 0; i < shortRead.length(); i++)
		{
			// mis match 발견 시
			if (refDNA[compare + i] != shortRead[i])
			{
				misMatchCount += 1;
			}
			// 발견된 mis match가 허용 수치를 넘길 경우
			if (misMatchCount > misMatchNum)
			{
				break;
			}
		}
		if (misMatchCount > misMatchNum)
		{
			return false;
		}
		else
		{
			// 비교 완료
			for (int i = 0; i < shortRead.length(); i++)
			{
				containerString[compare + i] = shortRead[i];
			}
			match = true;
			return true;
		}
	}
}

void findExactWithKMP::compareWithReference()
{
	ifstream mySequenceDNA("myDNA.txt");
	string comp;
	mySequenceDNA >> comp;
	int cnt01 = 0, cnt02 = 0;
	for (int i = 0; i < refDNA.length(); i++)
	{
		if (refDNA[i] != containerString[i])
			cnt01++;
		if (comp[i] != containerString[i])
			cnt02++;
	}

	cout << "reference DNA length : " << referenceDNASize << endl;
	cout << "reference mismatch : " << cnt01 << endl;
	cout << "my Sequence mismatch : " << cnt02 << endl;
	mySequenceDNA.close();
}


void findExactWithKMP::saveContainer()
{
	ofstream writeConstructed("constructed.txt");
	writeConstructed << containerString;
	writeConstructed.close();
}

int main()
{
	findExactWithKMP* aa = new findExactWithKMP();
	clock_t clockStart = clock();
	aa->run();
	aa->compareWithReference();
	double dueTime = ((float)clock() - clockStart) / CLOCKS_PER_SEC;
	int min = dueTime / 60;
	double sec = dueTime - 60 * min;
	printf("걸린 시간 :: %d분, %.2f 초 \n", min, sec);
	aa->saveContainer();
}

