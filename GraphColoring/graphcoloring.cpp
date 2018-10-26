// 禁忌搜索

#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
using namespace std;

const int MAX_COLOR_NUM = 50;  // 最大颜色数量
const int MAX_TEMP_SOLU = 50;  // 选择方案时预选的最大数量
const int MAX_RESTART_TIME = 10000; // 为找到可行方案时可以重开次数

struct Pair {
	int vertex;
	int color;
};

int vertexNum = 0;      // 顶点数量
int** adjacMatrix;  // 每个顶点的相邻边序号矩阵, 初始化时都为-1
int* adjacLen;  // 上一个矩阵每一行的实际大小
int* color;  // 顶点的颜色
int** funch;  // 对应颜色的顶点会有多少条边的两个顶点颜色相同
int** tabuIn;  // 禁忌表
Pair* increment;  // 邻域变换时不受禁忌点增量最大方案的数组
int increLen;  // 上一个数组的长度
Pair* incrementT;  // 邻域变换时受禁忌点增量最大方案的数组
int increTLen;  // 上一个数组的长度
int iter;  // 迭代次数
int f;  // 当前顶点颜色冲突数之和
int fBest; // 当前迭代最优f
int colorNum;

void readFile(string);  // 从文件中读取图的数据
bool tabu();
void outputForTest();
void output();
Pair selectMaxChange();
void update(Pair);  // 更新funch, color, f
void clearArray();  // 将数组funch, color, tabuIn的元素重置为0
int randomInt(int); //返回0~n的随机数
void initSolu();  // 随机初始解

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

//	 释放内存空间
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
		// 随机初始解
		initSolu();
		notImprove = 0;  // 记录迭代中解没有变优的次数

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
	int k = 0;  // 邻边序号

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
	int m = 0, m1 = -1, m2 = m1;  // m1为非禁忌方案最大增量，m2为禁忌方案
	for (int i = 0;i < vertexNum;i++) {
		if (funch[i][color[i]] != 0) {
			for (int j = 0;j < colorNum;j++) {
				// 不是当前颜色
				if (j != color[i]) {
					m = funch[i][color[i]] - funch[i][j];
					// 非禁忌方案
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
					} else {  // 禁忌方案
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

	// 特赦条件
	if (increTLen > 0 && m2 > m1 && (f - m2) < fBest) {
		m = randomInt(increTLen);
		return incrementT[m];
	} else if (increLen > 0) {  // 在非禁忌点中选择
		m = randomInt(increLen);
		return increment[m];
	} else {
		return {-1, -1};
	}
}

// 更新funch, color, f
void update(Pair p) {
	int k;  // 顶点序号
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

	//如果文件读取已经结束，表明当前文件为空文件，读取出错
	if (FIC.eof()) {
		cout << "### " << fileName << "is an empty file " << endl;
		exit(0);
	}

	int x1, x2; // 边对应的两个顶点
	int edgeNum, edgeIndex = 0;  // edgeIndex用于读取边时计数，判断边的数量与输入是否相符
	while (!FIC.eof()) {
		//如果StrReading第一个字符为p
		if (strcmp(StrReading, "p") == 0) {
			//从文件流中读取最大顶点数目和其边的数目，并输出
			FIC >> vertexNum >> edgeNum;
			cout << "Number of vertices = " << vertexNum << endl;
			cout << "Number of edges = " << edgeNum << endl;

			adjacMatrix = new int* [vertexNum];
			for (int i = 0;i < vertexNum;i++) {
				adjacMatrix[i] = new int [vertexNum];
				for (int j = 0;j < vertexNum;j++) {
					adjacMatrix[i][j] = -1;   // 顶点的编码可能为0，不能全部置为0
				}
			}
			adjacLen = new int[vertexNum] {0};
		}

		//如果StrReading第一个字符为e
		if (strcmp(StrReading, "e") == 0) {
			FIC >> x1 >> x2;
			x1--;x2--;
			// 检查输入顶点是否有误
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
	// 判断边数量与输入是否相符
	if (edgeIndex != edgeNum) {
		cout << "there is some problem with the egde num inputed: "
		     << "edgeIndex = " << edgeIndex << ", edgeNum = " << edgeNum << endl;
		exit(0);
	}
	cout << "Graph Density: " << edgeNum/(vertexNum*(vertexNum-1.0)) << endl;
	FIC.close();
	cout << "Finish reading data" << endl << endl;
}


