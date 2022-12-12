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
	string refDNA;	// ���۷��� DNA
	string shortReads; // ������
	string* shortReadsArray; // ������ �����带 ��� �迭
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

	cout << "�������� ���̿� ������ Ȯ��(%)�� �Է����ּ��� :: ";
	cin >> shortReadSize >> snpPercentage;
	dividedShortReadSize = shortReadSize / (snpPercentage + 1); // ������ �ɰ���
	misMatchNum = (shortReadSize / 100) * (snpPercentage) * 2; // ��� mismatch ��
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
	cout << endl << "KMP Time Count : " << milli01.count() << "milli second" << endl;
	cout << "�������� ���� : " << shortReadSize;
	cout << "���� �������� ���� : " << i << "��" << endl;
	cout << "��� mismatch �� : " << misMatchNum << endl;
}


void findExactWithKMP::divideAndCompare(string shortRead) {
	int position; // ������ �����尡 ���� �������� � ��ġ�� �ִ���

	shortReadsArray = new string[(snpPercentage + 1)]; // ������ ������ ���ڿ��� ���� �迭.
	match = false;

	// ���� �����带 ���� ��� �ݺ���
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		shortReadsArray[i] = shortRead.substr(i * dividedShortReadSize, dividedShortReadSize);
	}

	// ������ ��������� �ϳ��ϳ� �˻���. ����Ʈ ��Ī�� ���ϴ� ������.
	for (int i = 0; i < (snpPercentage + 1); i++)
	{
		SP.clear();
		// �ؽ�Ʈ�� ���̿� ������ ����
		int patternSize = dividedShortReadSize;
		int textSize = referenceDNASize;
		// ������ ���� ��ŭ �迭�� ����
		SP.assign(patternSize, 0);
		int j = 0;

		// i�� 1���� ����, j�� 0���� ����
		for (int k = 1; k < patternSize; k++)
		{
			// ���� j�� 0���� ũ�鼭, ������ i��ġ�� j��ġ�� �ٸ� ��� j�� �ڷ� �ϳ��� �̵�
			while (j > 0 && shortReadsArray[i][k] != shortReadsArray[i][j])
			{
				j = SP[j - 1];
			}

			// j�� 1����, i��ġ�� �߰�.
			if (shortReadsArray[i][k] == shortReadsArray[i][j])
			{
				SP[k] = ++j;
			}
		}

		// SP ���� �Ϸ�
		j = 0;

		int index = this->KMP(refDNA, shortReadsArray[i], 0, SP);

		if (index != -1)
		{
			while (!checkOtherStrings(index, i, shortRead))
			{
				index = this->KMP(refDNA, shortReadsArray[i], index + 1, SP);
				// �ٸ� �ε��� �˻�
				if (index == -1) break; // �˻� ����.
			}
			if (match)
			{
				// �ε����� �߰����� ��
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
	// 0 ���� �ؽ�Ʈ ���� ������
	for (int i = start; i < referenceDNASize; i++)
	{
		// ���λ� ���̻� ���̺��� j �ε����� ���ڰ� �ؽ�Ʈ�� i �ε����� ���ڿ� ��ġ�� �� ����.
		while (j > 0 && Text[i] != Pattern[j])
		{
			j = SP[j - 1];
		}

		// ���� ���ڰ� ���� ���
		if (Text[i] == Pattern[j])
		{
			// ���ڰ� ��� ��ġ�� ���
			if (j == dividedShortReadSize - 1)
			{
				return i - dividedShortReadSize + 1;
			}
			// ���ڰ� ��ġ�ϴ� ��� ���� ������ ����, j++
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
	else if (compare > 0) // �ش� ��ġ���� �����带 ���� �� �ִ���
	{
		for (int i = 0; i < shortRead.length(); i++)
		{
			// mis match �߰� ��
			if (refDNA[compare + i] != shortRead[i])
			{
				misMatchCount += 1;
			}
			// �߰ߵ� mis match�� ��� ��ġ�� �ѱ� ���
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
			// �� �Ϸ�
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
	printf("�ɸ� �ð� :: %d��, %.2f �� \n", min, sec);
	aa->saveContainer();
}

