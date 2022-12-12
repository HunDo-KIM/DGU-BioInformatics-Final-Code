#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <chrono>

// 2015112147 �赵��
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
	string refDNA;	// ���۷��� DNA
	string shortReads; // ������
	string* shortReadsArray; // ������ �����带 ��� �迭
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


// Bad Character ���̺� ����
void findExactWithBM::makeBadChar(string shortRead)
{
	// �ƽ�Ű �ڵ� ��ŭ ���� �Ҵ�
	BCContainer = new int[256];

	// ���̺� �� ��� ��ҵ��� �ʱ�ȭ
	for (int i = 0; i < 256; i++)
		BCContainer[i] = -1;

	// ���ڵ� ���� ���� ���. �տ� �ִ� ���ڰ� �ڿ��� �ߺ����� �߻� �� ���� ����
	for (int i = 0; i < dividedShortReadSize - 1; i++)
		BCContainer[(int)(shortRead[i])] = i;
}

// Good Suffix ���̺� ����
void findExactWithBM::makeGoodSuffix(string shortRead)
{
	int size = shortRead.length();
	GSContainer = new int[size];

	int i; // ���̺� ������ ��ġ
	int j = 0; 

	// ���� ���̺� �ʱ�ȭ, ������ ���̷� �ʱ�ȭ
	for (i = 0; i < size; i++)
	{
		GSContainer[i] = size;
	}

	i = size - 1; // ������ ũ�� M�� ������ �ε���, 0���� �����ϹǷ� M-1
	while (i > 0)
	{
		if (j >= 0 && shortRead[i + j] == shortRead[j])
			j -= 1;
		else 
		{
			if (j < 0)
			{
				// ���� ���̺ο� �´� ������ ã���� ���
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
	// i = �� ���� ����. j = ���� ��ġ
	int j = 0;
	int i = start;

	// ���̾��� �񱳸� �����ʺ���, i+j���� ��. ��ġ�ϸ� j ����
	while (i < (referenceDNASize - dividedShortReadSize))
	{
		j = dividedShortReadSize - 1;

		// �������� �������� �� ����
		while (j >= 0 && dividedShortRead[j] == refDNA[i + j])
			j--;
		if (j < 0) // j�� ������ ���� = ������ �߰��ߴٴ� ��
			return i;
		else // Good Suffix VS Bad Char, ���� ū������ ����Ʈ��.
			i += max(j - BCContainer[refDNA[i + j]], GSContainer[j]);
	}

	return -1;
}

void findExactWithBM::run()
{
	string inLine;
	int position, i = 0;
	cout << "�������� ���̿� ������ Ȯ��(���� %)�� �Է����ּ��� :: ";
	cin >> shortReadSize >> snpPercentage ;

	// �ΰ� DNA�� ���۷����� 3% ~4% ���̸� ����. 100�� ��� ���� ������ 
	dividedShortReadSize = shortReadSize / (snpPercentage + 1);
	misMatchNum = (shortReadSize / 100) * (snpPercentage)* 2;

	cout << "�м� ����" << endl;
	cout << "�м� ���� " << endl;
	chrono::system_clock::time_point StartTime01 = chrono::system_clock::now();
	while (getline(readShortReads, inLine))
	{
		cout << i << endl;
		divideAndCompare(inLine); // KMP ������ ��
		i++;
	}
	chrono::system_clock::time_point EndTime01 = chrono::system_clock::now();
	chrono::milliseconds milli01 = std::chrono::duration_cast<chrono::milliseconds>(EndTime01 - StartTime01);
	cout << endl << "BM Time Count : " << milli01.count() << "milli second" << endl;
	cout << "�������� ���� : " << shortReadSize;
	cout << "���� �������� ���� : " << i << "��" << endl;
	cout << "��� mismatch �� : " << misMatchNum << endl;
}

void findExactWithBM::divideAndCompare(string shortRead)
{
	int position; // ������ �����尡 ������ ��� �ִ���
	shortReadsArray = new string[(snpPercentage + 1)]; // ������ ������ ��Ʈ���� ���� �迭
	match = false;

	// ���� �����带 �����ϴ� �ݺ���
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		shortReadsArray[i] = shortRead.substr(i * dividedShortReadSize, dividedShortReadSize);
	}
	// ������ ��������� �ϳ��ϳ� �˻�. ����Ʈ ��Ī�� ���ϴ� ����
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
	printf("�ɸ� �ð� :: %d��, %.2f�� \n", min, sec);
	aa->saveContainer();

}