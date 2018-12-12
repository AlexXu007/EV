/*
	文件名：CppGATSP
	作者：Alex Xu
	地址：Dalian Maritime University
	描述：利用遗传算法求解TSP问题
	创建时间：2018年12月10日11点27分
*/

#include<iostream>
#include<vector>
#include<numeric>		//accumulate
#include<chrono>		//time
#include "HeuristicOperator.h"
using namespace std;
using namespace chrono;

//设置算法参数
# define	POP_SIZE	2
# define	MAX_GEN		4000

int main() {
	//计时开始
	auto start = system_clock::now();

	//生成距离矩阵
	HeuristicOperator ga_dm;
	vector<vector<double>> GA_DM;
	GA_DM = ga_dm.getDM(ga_dm.getCoord());

	int n = int(GA_DM[0].size());	//城市规模

	//初始化算法
	vector<vector<int>> initPop(POP_SIZE, vector<int>(n));		//初始种群
	vector<vector<int>> Pop(POP_SIZE, vector<int>(n));		//当前种群
	vector<vector<int>> newPop(POP_SIZE, vector<int>(n));		//新种群
	vector<double> popFit(POP_SIZE);		//记录种群适应度值
	vector<int> bestIndival(n);	//最优个体
	vector<double> gs(MAX_GEN + 1);	//记录全局最优解
	gs[0] = 1e9;
	unsigned int seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();

	//生成初始种群
	HeuristicOperator s0;
	for (int i = 0; i < POP_SIZE; i++) {
		initPop[i] = s0.getInitS(n);
	}
	Pop = initPop;

	//开始进化
	for (int gen = 1; gen <= MAX_GEN; gen++) {

		HeuristicOperator eval;			//计算种群的适应度值（这里直接用路径长度表示）
		for (int i = 0; i < POP_SIZE; i++) {
			popFit[i] = eval.Eval(Pop[i], GA_DM, n);
		}

		HeuristicOperator bestEI;		//找出种群中个体的最优适应度值并记录相应的个体编号
		vector<double> bestEvalIndex(2);
		bestEvalIndex = bestEI.bestS(popFit, POP_SIZE);
		double bestEval = bestEvalIndex[0];		//最优适应度值
		int bestIndex = int(bestEvalIndex[1]);	//最优适应度值对应的个体编号

		//最优解的更新
		if (bestEval < gs[gen - 1]) {		//比上一代优秀则更新
			gs[gen] = bestEval;
			bestIndival = Pop[bestIndex];
		}
		else {								//不比上一代优秀则不更新
			gs[gen] = gs[gen - 1];
		}
		if (gen % 1000 == 0) {
			cout << "第" << gen << "次迭代时全局最优评价值为" << gs[gen] << endl;
		}

		//扰动操作（产生新种群）
		for (int p = 0; p < POP_SIZE; p++) {
			HeuristicOperator shk;
			vector<int> randPosition = shk.RandPosition(n);
			vector<int> tmpS(n);
			double randShk = rand() / double(RAND_MAX);
			if (randShk < 0.33) {
				tmpS = shk.Swap(Pop[p], randPosition);		//交换操作
			}
			else if (randShk >= 0.67) {
				tmpS = shk.Flip(Pop[p], randPosition);		//翻转操作
			}
			else {
				tmpS = shk.Insert(Pop[p], randPosition);	//插入操作
			}

			HeuristicOperator evl;
			if (evl.Eval(tmpS, GA_DM, n) > evl.Eval(Pop[p], GA_DM, n)) {
				newPop[p] = Pop[p];
			}
			else {
				newPop[p] = tmpS;
			}
		}
		Pop = newPop;

		//选择操作（轮盘赌）
		vector<double> Cusum(POP_SIZE + 1, 0);		//适用于轮盘赌的累加器Cusum（补充了cus[0]=0;
		for (int i = 0; i < POP_SIZE; i++) {
			Cusum[i + 1] = Cusum[i] + popFit[i];
		}

		double Sum = accumulate(popFit.begin(), popFit.end(), 0.0);		//计算各个体被选择的概率（归一化）
		vector<double> cusFit(POP_SIZE + 1);		//放置种群中个个体被选择的概率值
		for (int i = 0; i < POP_SIZE + 1; i++) {
			cusFit[i] = Cusum[i] / Sum;
		}

		for (int p = 0; p < POP_SIZE; p++) {		//轮盘赌产生新种群
			double r = rand() / double(RAND_MAX);
			for (int q = 0; q < POP_SIZE; q++) {
				if (r > cusFit[q] && r <= cusFit[q + 1]) {
					newPop[p] = Pop[q];
				}
			}
		}
		Pop = newPop;
	}

	//计时结束
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << "花费了"
		<< double(duration.count()) * microseconds::period::num / microseconds::period::den
		<< "秒" << endl;

	//输出结果
	double gs0 = 15377.711;
	cout << "最优解为" << gs[MAX_GEN] << endl;
	double e = (gs[MAX_GEN] - gs0) / gs0;
	cout << "误差为" << e * 100.0 << '%' << endl;
	cout << "最优路径为" << endl;
	for (int i = 0; i < n; i++) {
		cout << bestIndival[i] + 1 << '\t';
	}
}