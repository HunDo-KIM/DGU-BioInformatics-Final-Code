#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <chrono>

// 2015112147 김도훈
using namespace std;

class findExactWithBM
{
public:
	findExactWithBM()
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
	~findExactWithBM()
	{
		readRefDNA.close();
		readShortReads.close();
	}

	int BM(string shortRead, int start);
	void makeGoodSuffix(string shortRead);
	void makeBadChar(string shortRead);

	void divideAndCompare(string shortRead);
	bool checkOtherStrings(int index, int position, string shortRead);
	void run();
	void compareWithRef();
	void saveContainer();


private:
	string refDNA;	// 레퍼런스 DNA
	string shortReads; // 숏리드
	string* shortReadsArray; // 나눠진 숏리드를 담는 배열
	string containerString;


	int shortReadSize;
	int dividedShortReadSize;
	int snpPercentage;
	int misMatchNum;
	int referenceDNASize;

	int* BCContainer;
	int* GSContainer;

	bool match;

	ifstream readRefDNA;
	ifstream readShortReads;
};


// Bad Character 테이블 생성
void findExactWithBM::makeBadChar(string shortRead)
{
	// 아스키 코드 만큼 동적 할당
	BCContainer = new int[256];

	// 테이블 내 모든 요소들을 초기화
	for (int i = 0; i < 256; i++)
		BCContainer[i] = -1;

	// 문자들 실제 출현 기록. 앞에 있던 문자가 뒤에서 중복으로 발생 시 덮어 씌움
	for (int i = 0; i < dividedShortReadSize - 1; i++)
		BCContainer[(int)(shortRead[i])] = i;
}

// Good Suffix 테이블 생성
void findExactWithBM::makeGoodSuffix(string shortRead)
{
	int size = shortRead.length();
	GSContainer = new int[size];

	int i; // 접미부 패턴의 위치
	int j = 0; 

	// 착한 접미부 초기화, 숏리드 길이로 초기화
	for (i = 0; i < size; i++)
	{
		GSContainer[i] = size;
	}

	i = size - 1; // 패턴의 크기 M의 마지막 인덱스, 0부터 시작하므로 M-1
	while (i > 0)
	{
		if (j >= 0 && shortRead[i + j] == shortRead[j])
			j -= 1;
		else 
		{
			if (j < 0)
			{
				// 착한 접미부와 맞는 패턴을 찾았을 경우
				while (i + j >= 0)
				{
					GSContainer[i + j] = i;
					j -= 1;
				}
			}
			else
				GSContainer[i + j] = i;
			j = size - i;
			i--;
		}
	}
}

int findExactWithBM::BM(string dividedShortRead, int start)
{
	// i = 비교 시작 지점. j = 패턴 위치
	int j = 0;
	int i = start;

	// 보이어무어는 비교를 오른쪽부터, i+j부터 비교. 일치하면 j 감소
	while (i < (referenceDNASize - dividedShortReadSize))
	{
		j = dividedShortReadSize - 1;

		// 우측에서 좌측으로 비교 시작
		while (j >= 0 && dividedShortRead[j] == refDNA[i + j])
			j--;
		if (j < 0) // j가 끝까지 갔다 = 패턴을 발견했다는 뜻
			return i;
		else // Good Suffix VS Bad Char, 둘중 큰것으로 쉬프트함.
			i += max(j - BCContainer[refDNA[i + j]], GSContainer[j]);
	}

	return -1;
}

void findExactWithBM::run()
{
	string inLine;
	int position, i = 0;
	cout << "숏리드의 길이와 스닙의 확률(정수 %)을 입력해주세요 :: ";
	cin >> shortReadSize >> snpPercentage ;

	// 인간 DNA가 레퍼런스와 3% ~4% 차이를 보임. 100의 배수 길이 숏리드 
	dividedShortReadSize = shortReadSize / (snpPercentage + 1);
	misMatchNum = (shortReadSize / 100) * (snpPercentage)* 2;

	cout << "분석 시작" << endl;
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
	cout << endl << "BM Time Count : " << milli01.count() << "milli second" << endl;
	cout << "숏리드의 길이 : " << shortReadSize;
	cout << "읽은 숏리드의 갯수 : " << i << "개" << endl;
	cout << "허용 mismatch 수 : " << misMatchNum << endl;
}

void findExactWithBM::divideAndCompare(string shortRead)
{
	int position; // 나눠진 숏리드가 원본의 어디에 있는지
	shortReadsArray = new string[(snpPercentage + 1)]; // 나눠진 숏리드 스트링을 담을 배열
	match = false;

	// 원본 숏리드를 분할하는 반복문
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		shortReadsArray[i] = shortRead.substr(i * dividedShortReadSize, dividedShortReadSize);
	}
	// 나눠진 숏리드들을 하나하나 검색. 퍼펙트 매칭을 구하는 과정
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		makeGoodSuffix(shortReadsArray[i]);
		makeBadChar(shortReadsArray[i]);
		int index = this->BM(shortReadsArray[i], 0);

		if (index != -1)
		{
			while (!checkOtherStrings(index, i, shortRead))
			{
				index = this->BM(shortReadsArray[i], index + 1);
				if (index == -1) break;
			}
			if (match)
				break;
			else
				continue;
		}
	}
	if (BCContainer != NULL && GSContainer != NULL)
	{
		delete[] shortReadsArray;
		delete[] BCContainer;
		delete[] GSContainer;
	}
}

bool findExactWithBM::checkOtherStrings(int index, int position, string shortRead)
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
		else
		{
			for (int i = 0; i < shortRead.length(); i++)
				containerString[compare + i] = shortRead[i];
			match = true;
			return true;
		}
	}
}

void findExactWithBM::saveContainer()
{
	ofstream writeConstructed("constructed.txt");
	writeConstructed << containerString;
	writeConstructed.close();
}

void findExactWithBM::compareWithRef()
{
	ifstream mySequenceDNA("myDNA.txt");
	string comp;
	mySequenceDNA >> comp;
	int cnt01 = 0, cnt02 = 0;
	for(int i = 0; i <refDNA.length(); i++)
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
	findExactWithBM* aa = new findExactWithBM();
	clock_t clockStart = clock();
	aa->run();
	aa->compareWithRef();
	double dueTime = ((float)clock() - clockStart) / CLOCKS_PER_SEC;
	int min = dueTime / 60;
	double sec = dueTime - 60 * min;
	printf("걸린 시간 :: %d분, %.2f초 \n", min, sec);
	aa->saveContainer();

}