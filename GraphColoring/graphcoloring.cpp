// ��������

#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
using namespace std;

const int MAX_COLOR_NUM = 50;  // �����ɫ����
const int MAX_TEMP_SOLU = 50;  // ѡ�񷽰�ʱԤѡ���������
const int MAX_RESTART_TIME = 10000; // Ϊ�ҵ����з���ʱ�����ؿ�����

struct Pair {
	int vertex;
	int color;
};

int vertexNum = 0;      // ��������
int** adjacMatrix;  // ÿ����������ڱ���ž���, ��ʼ��ʱ��Ϊ-1
int* adjacLen;  // ��һ������ÿһ�е�ʵ�ʴ�С
int* color;  // �������ɫ
int** funch;  // ��Ӧ��ɫ�Ķ�����ж������ߵ�����������ɫ��ͬ
int** tabuIn;  // ���ɱ�
Pair* increment;  // ����任ʱ���ܽ��ɵ�������󷽰�������
int increLen;  // ��һ������ĳ���
Pair* incrementT;  // ����任ʱ�ܽ��ɵ�������󷽰�������
int increTLen;  // ��һ������ĳ���
int iter;  // ��������
int f;  // ��ǰ������ɫ��ͻ��֮��
int fBest; // ��ǰ��������f
int colorNum;

void readFile(string);  // ���ļ��ж�ȡͼ������
bool tabu();
void outputForTest();
void output();
Pair selectMaxChange();
void update(Pair);  // ����funch, color, f
void clearArray();  // ������funch, color, tabuIn��Ԫ������Ϊ0
int randomInt(int); //����0~n�������
void initSolu();  // �����ʼ��

int main() {
	string fileName = "DSJC505.txt";
	readFile(fileName);

	color = new int[vertexNum];
	funch = new int* [vertexNum];
	tabuIn = new int* [vertexNum];
	for (int i = 0; i < vertexNum;i++) {
		funch[i] = new int [MAX_COLOR_NUM];
		tabuIn[i] = new int [MAX_COLOR_NUM];
	}
	increment = new Pair [MAX_TEMP_SOLU];
	incrementT = new Pair [MAX_TEMP_SOLU];
	colorNum = MAX_COLOR_NUM;

	while (colorNum > 1) {
		if (!tabu()) {
			cout << "colorNum = " << colorNum << endl;
			break;
		}
		colorNum--;
	}

//	 �ͷ��ڴ�ռ�
	for (int i = 0;i < vertexNum;i++) {
		delete[] adjacMatrix[i];
		delete[] funch[i];
		delete[] tabuIn[i];
	}
	delete[] adjacMatrix;
	delete[] adjacLen;
	delete[] funch;
	delete[] color;
	delete[] tabuIn;
	delete[] increment;
	delete[] incrementT;
	return 0;
}

bool tabu() {
	long startTime = clock();
	iter = 0;
	Pair k;
	int notImprove;

	for (int i = 0;i < MAX_RESTART_TIME;i++) {
		clearArray();
		// �����ʼ��
		initSolu();
		notImprove = 0;  // ��¼�����н�û�б��ŵĴ���

		while (f > 0 && notImprove < 100000) {
			k = selectMaxChange();

			if (k.vertex != -1) {
				update(k);
				if (fBest > f) {
					fBest = f;
					notImprove = 0;
				}
			} else {
				break;
			}
			iter++;
			notImprove++;
		}
		if (fBest == 0) {
			cout << endl;
			output();
			cout << "time used = " << double(clock() - startTime) / CLOCKS_PER_SEC << endl << endl;
			return true;
		} else {
			output();
			cout << "### restart" << endl;
		}
	}
	return false;
}

void initSolu() {
	int k = 0;  // �ڱ����

	f = 0;
	for (int i = 0;i < vertexNum;i++) {
		color[i] = randomInt(colorNum);
	}

	for (int i = 0;i < vertexNum;i++) {
		for (int j = 0;j < adjacLen[i];j++) {
			k = adjacMatrix[i][j];
			funch[i][color[k]]++;
		}
	}

	for (int i = 0;i < vertexNum;i++) {
		f += funch[i][color[i]];
	}
	f = f / 2;

	fBest = f;
}

Pair selectMaxChange() {
	increLen = increTLen = 0;
	int m = 0, m1 = -1, m2 = m1;  // m1Ϊ�ǽ��ɷ������������m2Ϊ���ɷ���
	for (int i = 0;i < vertexNum;i++) {
		if (funch[i][color[i]] != 0) {
			for (int j = 0;j < colorNum;j++) {
				// ���ǵ�ǰ��ɫ
				if (j != color[i]) {
					m = funch[i][color[i]] - funch[i][j];
					// �ǽ��ɷ���
					if (tabuIn[i][j] <= iter) {
						if (m > m1) {
							increLen = 0;
							m1 = m;
							increment[increLen] = {i, j};
							increLen++;
						} else if (m == m1 && increLen < MAX_TEMP_SOLU) {
							increment[increLen] = {i, j};
							increLen++;
						}
					} else {  // ���ɷ���
						if (m > m1 && m >= m2) {
							if (m > m2) {
								increTLen = 0;
								m2 = m;
								incrementT[increTLen] = {i, j};
								increTLen++;
							} else if (m == m2 && increTLen < MAX_TEMP_SOLU) {
								incrementT[increTLen] = {i, j};
								increTLen++;
							}
						}
					}
				}
			}
		}
	}

	// ��������
	if (increTLen > 0 && m2 > m1 && (f - m2) < fBest) {
		m = randomInt(increTLen);
		return incrementT[m];
	} else if (increLen > 0) {  // �ڷǽ��ɵ���ѡ��
		m = randomInt(increLen);
		return increment[m];
	} else {
		return {-1, -1};
	}
}

// ����funch, color, f
void update(Pair p) {
	int k;  // �������
	for (int i = 0;i < adjacLen[p.vertex];i++) {
		k = adjacMatrix[p.vertex][i];
		funch[k][p.color]++;
		funch[k][color[p.vertex]]--;
	}
	f -= funch[p.vertex][color[p.vertex]] - funch[p.vertex][p.color];
	tabuIn[p.vertex][color[p.vertex]] = iter + f + randomInt(10);
	color[p.vertex] = p.color;
}

void outputForTest() {
	cout << "iter time: " << iter << "; f = " << f <<endl;

//	cout << endl << "adjaMatrix:  " << endl;
//	for (int i = 0;i < vertexNum;i++) {
//		for (int j = 0;j < vertexNum;j++) {
//			cout << adjacMatrix[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << endl << "adjacLen: " << endl;
//	for (int i = 0;i < vertexNum;i++) {
//		cout << adjacLen[i] << " ";
//	}
//	cout << endl;

	cout << endl << "color: " << endl;
	for (int i = 0;i < vertexNum;i++) {
		cout << color[i] << " ";
	}
	cout << endl;

	cout << endl << "funch:  " << endl;
	for (int i = 0;i < vertexNum;i++) {
		for (int j = 0;j < colorNum;j++) {
			cout << funch[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl << "tabuIn:  " << endl;
	for (int i = 0;i < vertexNum;i++) {
		for (int j = 0;j < colorNum;j++) {
			cout << tabuIn[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl << "increment: " << endl;
	for (int i = 0;i < increLen;i++) {
		cout << increment[i].vertex << "," << increment[i].color << "  ";
	}
	cout << endl << "incrementT: " << endl;
	for (int i = 0;i < increTLen;i++) {
		cout << incrementT[i].vertex << "," << incrementT[i].color << "  ";
	}
	cout << endl;
}

void output() {
	cout << "colorNum = " << colorNum << "; fBest = " << fBest << "; iter time = " << iter << endl;
}

void clearArray() {
	for (int i = 0;i < vertexNum;i++) {
		memset(funch[i], 0, sizeof(int) * colorNum);
		memset(tabuIn[i], 0, sizeof(int) * colorNum);
	}
}

int randomInt(int n) {
	return rand() % n;
}

void readFile(string fileName) {
	ifstream FIC;

	FIC.open(fileName);
	if (FIC.fail()) {
		cout << "### Error open file " << fileName << endl;
		getchar();
		exit(0);
	}

	char StrReading[100];
	FIC >> StrReading;

	//����ļ���ȡ�Ѿ�������������ǰ�ļ�Ϊ���ļ�����ȡ����
	if (FIC.eof()) {
		cout << "### " << fileName << "is an empty file " << endl;
		exit(0);
	}

	int x1, x2; // �߶�Ӧ����������
	int edgeNum, edgeIndex = 0;  // edgeIndex���ڶ�ȡ��ʱ�������жϱߵ������������Ƿ����
	while (!FIC.eof()) {
		//���StrReading��һ���ַ�Ϊp
		if (strcmp(StrReading, "p") == 0) {
			//���ļ����ж�ȡ��󶥵���Ŀ����ߵ���Ŀ�������
			FIC >> vertexNum >> edgeNum;
			cout << "Number of vertices = " << vertexNum << endl;
			cout << "Number of edges = " << edgeNum << endl;

			adjacMatrix = new int* [vertexNum];
			for (int i = 0;i < vertexNum;i++) {
				adjacMatrix[i] = new int [vertexNum];
				for (int j = 0;j < vertexNum;j++) {
					adjacMatrix[i][j] = -1;   // ����ı������Ϊ0������ȫ����Ϊ0
				}
			}
			adjacLen = new int[vertexNum] {0};
		}

		//���StrReading��һ���ַ�Ϊe
		if (strcmp(StrReading, "e") == 0) {
			FIC >> x1 >> x2;
			x1--;x2--;
			// ������붥���Ƿ�����
			if (x1<0 || x2<0 || x1>=vertexNum || x2 >=vertexNum ) {
				cout << "### Error of edge code : x1="
						 << x1 << ", x2=" << x2 << endl;
				exit(0);
			}
			adjacMatrix[x1][adjacLen[x1]] = x2; adjacLen[x1]++;
			adjacMatrix[x2][adjacLen[x2]] = x1; adjacLen[x2]++;
			edgeIndex++;
		}
		FIC >> StrReading;
	}
	// �жϱ������������Ƿ����
	if (edgeIndex != edgeNum) {
		cout << "there is some problem with the egde num inputed: "
		     << "edgeIndex = " << edgeIndex << ", edgeNum = " << edgeNum << endl;
		exit(0);
	}
	cout << "Graph Density: " << edgeNum/(vertexNum*(vertexNum-1.0)) << endl;
	FIC.close();
	cout << "Finish reading data" << endl << endl;
}


