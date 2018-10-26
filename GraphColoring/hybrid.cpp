#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

const int P_SIZE = 10;  // 种群大小
int iterL = 5600;  // 禁忌搜索迭代次数
const int MAX_COLOR_NUM = 48;
const int MAX_TEMP_SOLU = 50;  // 禁忌搜索中选择方案时预选的最大数量
const int MAX_F = 9999999;
const int MAX_RESTART = 100;   // 最大重启次数

struct Pair {
	int vertex;
	int color;
};

int vertexNum;  // 顶点数目
int colorNum;  // 所用颜色数目

int** graphMatrix;  // 图的邻接矩阵
int** adjacMatrix;  // 保存顶点的相邻边
int* lenAdjac;  // 上一个矩阵每一行的实际大小

int** population;  // 种群，其中一个保存子代
int** lenPopu;  // 种群中个体每个子集的大小，其中一个保存子代
int** address;  // 种群中个体的顶点位置，其中一个保存子代的顶点位置,-1为未分配，-2为分配失败
int* fPopu;  // 种群中个体的有冲突的边数
int childIndex;  // 上面三个数组中子代的序号

int* temp;    // 标志交叉互换时是否已选
int* tempP1;  // 保存交叉互换时P1各集合的长度
int* tempP2;  // 保存交叉互换时P1各集合的长度
int* tempM;   //  保存颜色集中数量最大的集合序号

int* color;
int* C_color;  //当前最优解
int** funch;       // 对应颜色的顶点会有多少条边的两个顶点颜色相同
int** tabuIn;      // 禁忌表
Pair* increment;   // 邻域变换时不受禁忌点增量最大方案的数组
int increLen;      // 上一个数组的长度
Pair* incrementT;  // 邻域变换时受禁忌点增量最大方案的数组
int increTLen;     // 上一个数组的长度
int iter;
int fBest;
int f;
int crossTime;   // 交叉互换次数
Pair minF;

void initPopulation();  // 初始化种群
void crossover(int, int);  //  交叉互换，传入参数为父代在population中的序号
void localSearch(int);  // 局部搜索优化解
Pair selectMaxChange();
void update(Pair);
void initSolu(int);
void clearArray();   // reset funch and tabuIn
void resetArray();   // reset lenPopu, fPopu and address
void updatePopu(int, int);  // 更新种群
int hasSame(int, int);  // 判断种群是否有相同个体
int randomInt(int);
int findMax(int*, int);  // 返回最大值的序号
void randomInsert(int, int, int, int);
void readFile(string);
void outputForTest();
void output();

int main() {
	int p1, p2;
	int notImprove;
	colorNum = MAX_COLOR_NUM;

	readFile("DSJC505.txt");

	population = new int*[P_SIZE+1];
	lenPopu = new int*[P_SIZE+1];
	address = new int*[P_SIZE+1];
	for(int i = 0;i < P_SIZE+1;i++) {
		population[i] = new int[vertexNum];
		lenPopu[i] = new int[MAX_COLOR_NUM] {0};
		address[i] = new int[vertexNum];
		for (int j = 0;j < vertexNum;j++) {
			address[i][j] = -1;  // 未分配
		}
	}
	fPopu = new int[P_SIZE+1] {0};
	temp = new int[vertexNum];
	tempP1 = new int[colorNum];
	tempP2 = new int[colorNum];
	tempM = new int[colorNum];

	funch = new int* [vertexNum];
	tabuIn = new int* [vertexNum];
	for (int i = 0; i < vertexNum;i++) {
		funch[i] = new int [MAX_COLOR_NUM];
		tabuIn[i] = new int [MAX_COLOR_NUM];
	}
	increment = new Pair [MAX_TEMP_SOLU];
	incrementT = new Pair [MAX_TEMP_SOLU];
	color = new int[vertexNum];
	C_color = new int[vertexNum];

	for (int i = 0;i < MAX_RESTART;i++) {
		cout << "****  run time : " << i << endl << endl;
		resetArray();
		childIndex = P_SIZE;

		initPopulation();
		iter = 0;
		notImprove = 0;
		crossTime = 0;
		while (minF.vertex > 0 && notImprove < 10000 && crossTime < 2000) {
			// 随机选择父代
			p1 = randomInt(P_SIZE+1);
			if (p1 == childIndex) {
				p1 = (p1 + 1) % (P_SIZE+1);
			}
			// 保证两次随机数不会相同，且第二次随机不是子代
			do {
				p2 = (p1 + randomInt(P_SIZE) + 1) % (P_SIZE+1);
			} while (p2 == childIndex);
			// 交叉互换
			crossover(p1, p2);
			crossTime++;
			localSearch(childIndex);

			if (minF.vertex > fPopu[childIndex]) {
				minF = {fPopu[childIndex], childIndex};
				notImprove = 0;
				output();
			} else {
				notImprove++;
			}

			updatePopu(p1, p2);
		}
		if (minF.vertex == 0) {
			break;
		}
	}


	cout << endl << "**** final answer :" << endl;
	output();

	// 释放数组空间
	for (int i = 0;i < vertexNum;i++) {
		delete[] graphMatrix[i];
		delete[] adjacMatrix[i];
		delete[] funch[i];
		delete[] tabuIn[i];
	}
	delete[] graphMatrix;
	delete[] adjacMatrix;
	delete[] lenAdjac;
	delete[] funch;
	delete[] tabuIn;
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

	delete[] increment;
	delete[] incrementT;
	delete[] color;
	delete[] C_color;
	return 0;
}
 // reset lenPopu and fPopu and address
void resetArray() {
	for (int i = 0;i < P_SIZE+1;i++) {
		fPopu[i] = 0;
		for (int j = 0;j < colorNum;j++) {
			lenPopu[i][j] = 0;
		}
		for (int j = 0;j < vertexNum;j++) {
			address[i][j] = -1;
		}
	}
}

// 初始化种群
void initPopulation() {
	minF = {MAX_F, -1};
	int v, vi, vm, va, ca, n;  // va 为已指定的顶点个数, ca为已使用的颜色个数
	bool flag;
	for (int i = 0;i < P_SIZE;i++) {
		va = ca = 0;
		for (int j = 0;j < vertexNum;j++) {
			v = randomInt(vertexNum);
			while (address[i][v] != -1) {
				v = (v + 1) % vertexNum;
			}

			vi = 0;
			// 指定顶点的集合
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
					// 判断有无邻接点
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
					// 分配不成功的顶点
					if (k == colorNum-1 && flag == false) {
						address[i][v] = -2;
					}
				}
				vi += lenPopu[i][k];
			}
		}

		// 随机指定剩余的顶点
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

		// 优化解
		localSearch(i);

		if (hasSame(i, 30) != -1) {
			for (int j = 0;j < vertexNum;j++) {
				address[i][j] = -1;
			}
			for (int j = 0;j < colorNum;j++) {
				lenPopu[i][j] = 0;
			}
			i--;
		} else {
			if (minF.vertex > fPopu[i]) {
				minF = {fPopu[i], i};
			}
		}
	}
}

// 交叉互换，传入参数为父代在population中的序号
void crossover(int p1, int p2) {
	int t1, cLen = 0, si, v, t2, k, t3;
	for (int i = 0;i < vertexNum;i++) {
		temp[i] = 0;  // 标记顶点是否已选
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
				// 将选出的点复制到子代中
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
				// 将选出的点复制到子代中
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

	// 随机分配剩余点
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

// vi为插入的颜色序号，vm 为已有的顶点个数，i为插入个体在种群中的序号, v为顶点序号
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

// 更新种群
void updatePopu(int p1, int p2) {
	int vi = 0;
	// 可以得到比当前最优解更好的解
	if (fPopu[childIndex] < minF.vertex) {
		if (fPopu[p1] <= fPopu[p2]) {
			childIndex = p2;
		} else {
			childIndex = p1;
		}
	} else {
		vi = hasSame(childIndex, 6);
		// 有类似解但是比类似解更好，替换该解
		if (vi != -1) {
			if (fPopu[childIndex] < fPopu[vi]) {
				childIndex = vi;
			}
		} else {
			// 没有类似解判断解的质量决定是否替换
			if (fPopu[p1] < fPopu[p2]) {
				if (fPopu[childIndex] <= fPopu[p2]) {
					childIndex = p2;
				}
			} else if (fPopu[p1] == fPopu[p2]) {
				if (fPopu[childIndex] <= fPopu[p2]) {
					vi = randomInt(2);
					if (vi == 0)
						childIndex = p1;
					else
						childIndex = p2;
				}
			} else {
				if (fPopu[childIndex] <= fPopu[p1]) {
					childIndex = p1;
				}
			}
		}
	}
}

// 判断种群是否有类似个体
int hasSame(int pi, int bound) {
	int re = 0, c = 0;
	for (int i = 0;i < P_SIZE+1;i++) {
		if (address[i][0] != -1 && i != pi) {
			re = 0;
			for (int j = 0;j < vertexNum;j++) {
				if (population[i][j] != population[pi][j]) {
					re++;
				}
				if (re >= bound) {
					break;
				}
			}
			// 只有小于 bound 个数据不一样，且分组相同
			if (re < bound) {
				c = 0;
				for (int k = 0;k < colorNum;k++) {
					if (lenPopu[i][k] == lenPopu[pi][k]) {
						c++;
					}
				}
				if (c == colorNum) {
					return i;
				}
			}
		}
	}
	return -1;
}

// 局部搜索优化解,pi 为个体在种群中的序号
void localSearch(int pi) {
	int vi, ci;
	clearArray();

	initSolu(pi);
	Pair k;
	for (int i = 0;i < iterL;i++) {
		k = selectMaxChange();
		if (k.vertex != -1) {
			update(k);
			if (f <= fBest) {
				fBest = f;
				for (int j = 0;j < vertexNum;j++) {
					C_color[j] =  color[j];
				}
				if (fBest == 0) {
					break;
				}
			}
		} else {
			break;
		}
		iter++;
	}

	vi = 0;ci = 0;
	for (int i = 0;i < colorNum;i++) {
		for (int j = 0;j < vertexNum;j++) {
			if (C_color[j] == i) {
				population[pi][vi] = j;
				address[pi][j] = vi;
				vi++;
				ci++;
			}
		}
		lenPopu[pi][i] = ci;
		ci = 0;
	}
	fPopu[pi] = fBest;
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
	for (int i = 0;i < lenAdjac[p.vertex];i++) {
		k = adjacMatrix[p.vertex][i];
		funch[k][p.color]++;
		funch[k][color[p.vertex]]--;
	}
	f -= funch[p.vertex][color[p.vertex]] - funch[p.vertex][p.color];
	tabuIn[p.vertex][color[p.vertex]] = iter + int(0.6 * f) + randomInt(10);
	color[p.vertex] = p.color;
}

void initSolu(int pi) {
	int k = 0, vi = 0;  // 邻边序号

	f = 0;
	for (int i = 0;i < colorNum;i++) {
		for (int j = 0;j < lenPopu[pi][i];j++) {
			color[population[pi][vi]] = i;
			C_color[population[pi][vi]] = i;
			vi++;
		}
	}

	for (int i = 0;i < vertexNum;i++) {
		for (int j = 0;j < lenAdjac[i];j++) {
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

void clearArray() {
	for (int i = 0;i < vertexNum;i++) {
		memset(funch[i], 0, sizeof(int) * colorNum);
		memset(tabuIn[i], 0, sizeof(int) * colorNum);
	}
}

int randomInt(int n) {
	return rand() % n;
}

// 返回最大值的序号
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
			graphMatrix = new int* [vertexNum];
			for (int i = 0;i < vertexNum;i++) {
				adjacMatrix[i] = new int [vertexNum];
				for (int j = 0;j < vertexNum;j++) {
					adjacMatrix[i][j] = -1;   // 顶点的编码可能为0，不能全部置为0
				}
				graphMatrix[i] = new int[vertexNum] {0};
			}
			lenAdjac = new int[vertexNum] {0};
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
			adjacMatrix[x1][lenAdjac[x1]] = x2; lenAdjac[x1]++;
			adjacMatrix[x2][lenAdjac[x2]] = x1; lenAdjac[x2]++;

			graphMatrix[x1][x2] = 1;
			graphMatrix[x2][x1] = 1;
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

	cout << "population : " << endl;
	for (int i = 0;i < P_SIZE+1;i++) {
		for (int j = 0;j < vertexNum;j++) {
			cout << population[i][j] << " ";
		}
		cout << endl;
	}

	cout << "lenPopu : " << endl;
	for (int i = 0;i < P_SIZE+1;i++) {
		for (int j = 0;j < colorNum;j++) {
			cout << lenPopu[i][j] << " ";
		}
		cout << endl;
	}

	cout << "address : " << endl;
	for (int i = 0;i < P_SIZE+1;i++) {
		for (int j = 0;j < vertexNum;j++) {
			cout << address[i][j] << " ";
		}
		cout << endl;
	}

	cout << "fPopu : " << endl;
	for (int i = 0;i < P_SIZE+1;i++) {
		cout << fPopu[i] <<  " ";
	}
	cout << endl;

//	cout << "temp : " << endl;
//	for (int i = 0;i < vertexNum;i++) {
//		cout << temp[i] <<  " ";
//	}
//	cout << endl;
//
//	cout << "tempP1 : " << endl;
//	for (int i = 0;i < colorNum;i++) {
//		cout << tempP1[i] <<  " ";
//	}
//	cout << endl;
//
//	cout << "tempP2 : " << endl;
//	for (int i = 0;i < colorNum;i++) {
//		cout << tempP2[i] <<  " ";
//	}
//	cout << endl;

//	cout << "child : " << endl;
//	cout << "vertex: ";
//	for (int i = 0;i < vertexNum;i++) {
//		cout << population[childIndex][i] <<  " ";
//	}
//	cout << endl << "address: ";
//	for (int i = 0;i < vertexNum;i++) {
//		cout << address[childIndex][i] <<  " ";
//	}
//	cout << endl << "len: ";
//	for (int i = 0;i < colorNum;i++) {
//		cout << lenPopu[childIndex][i] <<  " ";
//	}
//	cout << endl;

	cout << endl << endl;
}

void output() {
	cout << endl;
	cout << "iter = " << iter << " ; crossTime = " << crossTime;
	cout << " ; f = " << minF.vertex ;
	cout << endl;
}
