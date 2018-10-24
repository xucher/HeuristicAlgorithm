#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

const int P_SIZE = 5;  // ��Ⱥ��С
const int MAX_COLOR_NUM = 3;

int vertexNum;  // ������Ŀ
int colorNum;  // ������ɫ��Ŀ

int** graphMatrix;  // ͼ���ڽӾ���
int** adjacMatrix;  // ���涥������ڱ�
int* lenAdjac;  // ��һ������ÿһ�е�ʵ�ʴ�С

int** population;  // ��Ⱥ������һ�������Ӵ�
int** lenPopu;  // ��Ⱥ�и���ÿ���Ӽ��Ĵ�С������һ�������Ӵ�
int** address;  // ��Ⱥ�и���Ķ���λ�ã�����һ�������Ӵ��Ķ���λ��,-1Ϊδ���䣬-2Ϊ����ʧ��
int* fPopu;  // ��Ⱥ�и�����г�ͻ�ı���
int childIndex;  // ���������������Ӵ������

int* temp;    // ��־���滥��ʱ�Ƿ���ѡ
int* tempP1;  // ���潻�滥��ʱP1�����ϵĳ���
int* tempP2;  // ���潻�滥��ʱP1�����ϵĳ���
int* tempM;   //  ������ɫ�����������ļ�����ţ�

void initPopulation();  // ��ʼ����Ⱥ
void crossover(int, int);  //  ���滥�����������Ϊ������population�е����
void localSearch();  // �ֲ������Ż���
void updatePopu();  // ������Ⱥ
int randomInt(int);
int findMax(int*, int);  // �������ֵ�����
void randomInsert(int, int, int, int);
void readFile(string);
void outputForTest();

int main() {
	int p1, p2;
	colorNum = MAX_COLOR_NUM;

	readFile("test.txt");

	population = new int*[P_SIZE+1];
	lenPopu = new int*[P_SIZE+1];
	address = new int*[P_SIZE+1];
	for(int i = 0;i < P_SIZE+1;i++) {
		population[i] = new int[vertexNum];
		lenPopu[i] = new int[MAX_COLOR_NUM] {0};
		address[i] = new int[vertexNum];
		for (int j = 0;j < vertexNum;j++) {
			address[i][j] = -1;  // δ����
		}
	}
	fPopu = new int[P_SIZE+1] {0};
	temp = new int[vertexNum];
	tempP1 = new int[colorNum];
	tempP2 = new int[colorNum];
	tempM = new int[colorNum];

	childIndex = P_SIZE;
	initPopulation();

	while (1) {
		// ���ѡ�񸸴�
		p1 = randomInt(P_SIZE+1);
		if (p1 == childIndex) {
			p1 = (p1 + 1) % (P_SIZE+1);
		}
		// ��֤���������������ͬ���ҵڶ�����������Ӵ�
		do {
			p2 = (p1 + randomInt(P_SIZE) + 1) % (P_SIZE+1);
		} while (p2 == childIndex);
		// ���滥��
		crossover(p1, p2);
		localSearch();
		updatePopu();
//		break;
	}
	// �ͷ�����ռ�
	for (int i = 0;i < vertexNum;i++) {
		delete[] graphMatrix[i];
		delete[] adjacMatrix[i];
	}
	delete[] graphMatrix;
	delete[] adjacMatrix;
	delete[] lenAdjac;
	for (int i = 0;i < P_SIZE;i++) {
		delete[] population[i];
		delete[] lenPopu[i];
		delete[] address[i];
	}
	delete[] population;
	delete[] lenPopu;
	delete[] address;
	delete[] fPopu;
	delete[] temp;
	delete[] tempP1;
	delete[] tempP2;
	delete[] tempM;
	return 0;
}

// ��ʼ����Ⱥ
void initPopulation() {
	int v, vi, vm, va, ca, n;  // va Ϊ��ָ���Ķ������, caΪ��ʹ�õ���ɫ����
	bool flag;
	for (int i = 0;i < P_SIZE;i++) {
		va = ca = 0;
		for (int j = 0;j < vertexNum;j++) {
			v = randomInt(vertexNum);
			while (address[i][v] != -1) {
				v = (v + 1) % vertexNum;
			}

			vi = 0;
			// ָ������ļ���
			for (int k = 0;k < colorNum;k++) {
				if (lenPopu[i][k] == 0) {
					population[i][vi] = v;
					lenPopu[i][k]++;
					address[i][v] = vi;
					va++;
					ca++;
					break;
				} else {
					flag = true;
					// �ж������ڽӵ�
					for (int l = 0;l < lenPopu[i][k];l++) {
						if (graphMatrix[v][population[i][vi+l]] == 1) {
							flag = false;break;
						}
					}
					if (flag == true) {
						vm = va; n = 1;
						while (ca - n > k) {
							population[i][vm] = population[i][vm-lenPopu[i][ca-n]];
							address[i][population[i][vm]] = vm;
							vm -= lenPopu[i][ca-n];
							n++;
						}
						population[i][vm] = v;
						lenPopu[i][k]++;
						address[i][v] = vm;
						va++;
						break;
					}
					// ���䲻�ɹ��Ķ���
					if (k == colorNum-1 && flag == false) {
						address[i][v] = -2;
					}
				}
				vi += lenPopu[i][k];
			}
		}

		// ���ָ��ʣ��Ķ���
		for (int j = 0;j < vertexNum;j++) {
			if (address[i][j] == -2) {
				vi = randomInt(colorNum);
				randomInsert(vi, va, i, j);
				va++;
				vm = 0;
				for (int k = 0;k < vi;k++) {
					vm += lenPopu[i][k];
				}

				for (int k = 0;k < lenPopu[i][vi]-1;k++) {
					if (graphMatrix[j][population[i][vm+k]] == 1) {
						fPopu[i]++;
					}
				}
			}
			if (va == vertexNum) {
				break;
			}
		}
	}
}

// ���滥�����������Ϊ������population�е����
void crossover(int p1, int p2) {
	int t1, cLen = 0, si, v, t2, k, t3;
	for (int i = 0;i < vertexNum;i++) {
		temp[i] = 0;  // ��Ƕ����Ƿ���ѡ
	}
	for (int i = 0;i < colorNum;i++) {
		tempP1[i] = lenPopu[p1][i];
		tempP2[i] = lenPopu[p2][i];
	}

	for (int i = 0;i < colorNum;i++) {
		if (i % 2 == 0) {
			t1 = findMax(tempP1, colorNum);
			si = 0;
			for (int j = 0;j < t1;j++) {
				si += lenPopu[p1][j];
			}

			for (int j = 0;j < lenPopu[p1][t1];j++) {
				// ��ѡ���ĵ㸴�Ƶ��Ӵ���
				v = population[p1][si+j];
				if (temp[v] == 0) {
					population[childIndex][cLen] = v;
					address[childIndex][v] = cLen;
					cLen++;
					temp[v] = 1;

					t2 = address[p2][v]; k = 0;t3 = 0;
					while (k < colorNum) {
						if ( t2 >= t3 && t2 < t3+lenPopu[p2][k]) {
							tempP2[k]--;
							break;
						}
						t3 += lenPopu[p2][k];
						k++;
					}
				}
			}
			lenPopu[childIndex][i] = tempP1[t1];
			tempP1[t1] = 0;
		} else {
			t1 = findMax(tempP2, colorNum);
			si = 0;
			for (int j = 0;j < t1;j++) {
				si += lenPopu[p2][j];
			}

			for (int j = 0;j < lenPopu[p2][t1];j++) {
				// ��ѡ���ĵ㸴�Ƶ��Ӵ���
				v = population[p2][si+j];
				if (temp[v] == 0) {
					population[childIndex][cLen] = v;
					address[childIndex][v] = cLen;
					cLen++;
					temp[v] = 1;

					t2 = address[p1][v]; k = 0;t3 = 0;
					while (k < colorNum) {
						if ( t2 >= t3 && t2 < t3+lenPopu[p1][k]) {
							tempP1[k]--;
							break;
						}
						t3 += lenPopu[p1][k];
						k++;
					}
				}
			}
			lenPopu[childIndex][i] = tempP2[t1];
			tempP2[t1] = 0;
		}
	}

	// �������ʣ���
  for (int i = 0;i < vertexNum;i++) {
		if (temp[i] == 0) {
			randomInsert(randomInt(colorNum), cLen, childIndex, i);
			cLen++;
		}

		if (cLen == vertexNum) {
			break;
		}
  }
}

// viΪ�������ɫ��ţ�vm Ϊ���еĶ��������iΪ�����������Ⱥ�е����, vΪ�������
void randomInsert(int vi, int vm, int i, int v) {
	int n = 1;
	while (colorNum - n > vi) {
		population[i][vm] = population[i][vm-lenPopu[i][colorNum-n]];
		address[i][population[i][vm]] = vm;
		vm -= lenPopu[i][colorNum-n];
		n++;
	}
	population[i][vm] = v;
	lenPopu[i][vi]++;
	address[i][v] = vm;
}

// �ֲ������Ż���
void localSearch() {

}

// ������Ⱥ
void updatePopu() {

}

int randomInt(int n) {
	return rand() % n;
}

// �������ֵ�����
int findMax(int* arr, int len) {
	int m = -1, mLen = 0;
	for (int i = 0;i < len;i++) {
		if (arr[i] > m) {
			m = arr[i];
			mLen = 0;
			tempM[mLen] = i;
			mLen++;
		} else if (arr[i] == m) {
			tempM[mLen] = i;
			mLen++;
		}
	}
	return tempM[randomInt(mLen)];
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
			graphMatrix = new int* [vertexNum];
			for (int i = 0;i < vertexNum;i++) {
				adjacMatrix[i] = new int [vertexNum];
				for (int j = 0;j < vertexNum;j++) {
					adjacMatrix[i][j] = -1;   // ����ı������Ϊ0������ȫ����Ϊ0
				}
				graphMatrix[i] = new int[vertexNum] {0};
			}
			lenAdjac = new int[vertexNum] {0};
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
			adjacMatrix[x1][lenAdjac[x1]] = x2; lenAdjac[x1]++;
			adjacMatrix[x2][lenAdjac[x2]] = x1; lenAdjac[x2]++;

			graphMatrix[x1][x2] = 1;
			graphMatrix[x2][x1] = 1;
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

void outputForTest() {
//	cout << "adjacMatrix : " << endl;
//	for (int i = 0;i < vertexNum;i++) {
//		for (int j = 0;j < lenAdjac[i];j++) {
//			cout << adjacMatrix[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << "graphMatrix : " << endl;
//	for (int i = 0;i < vertexNum;i++) {
//		for (int j = 0;j < vertexNum;j++) {
//			cout << graphMatrix[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << "population : " << endl;
//	for (int i = 0;i < P_SIZE+1;i++) {
//		for (int j = 0;j < vertexNum;j++) {
//			cout << population[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << "lenPopu : " << endl;
//	for (int i = 0;i < P_SIZE+1;i++) {
//		for (int j = 0;j < colorNum;j++) {
//			cout << lenPopu[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << "address : " << endl;
//	for (int i = 0;i < P_SIZE+1;i++) {
//		for (int j = 0;j < vertexNum;j++) {
//			cout << address[i][j] << " ";
//		}
//		cout << endl;
//	}
//
//	cout << "fPopu : " << endl;
//	for (int i = 0;i < P_SIZE+1;i++) {
//		cout << fPopu[i] <<  " ";
//	}
//	cout << endl;

	cout << "temp : " << endl;
	for (int i = 0;i < vertexNum;i++) {
		cout << temp[i] <<  " ";
	}
	cout << endl;

	cout << "tempP1 : " << endl;
	for (int i = 0;i < colorNum;i++) {
		cout << tempP1[i] <<  " ";
	}
	cout << endl;

	cout << "tempP2 : " << endl;
	for (int i = 0;i < colorNum;i++) {
		cout << tempP2[i] <<  " ";
	}
	cout << endl;

	cout << "child : " << endl;
	cout << "vertex: ";
	for (int i = 0;i < vertexNum;i++) {
		cout << population[childIndex][i] <<  " ";
	}
	cout << endl << "address: ";
	for (int i = 0;i < vertexNum;i++) {
		cout << address[childIndex][i] <<  " ";
	}
	cout << endl << "len: ";
	for (int i = 0;i < colorNum;i++) {
		cout << lenPopu[childIndex][i] <<  " ";
	}
	cout << endl;

	cout << endl << endl;
}

