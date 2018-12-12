#pragma once
#include<iostream>
#include<vector>
#include <algorithm>		// std::shuffle
#include <random>			// std::default_random_engine
#include<chrono>
using namespace std;

class HeuristicOperator {
public:
	vector<vector<double>> getCoord(void);		//函数1：获取坐标函数
	vector<vector<double>> getDM(vector<vector<double>> Coord);		//函数2：获取距离矩阵函数
	vector<int> getInitS(int n);		//函数3：获取初始解函数
	double Eval(vector<int> S, vector<vector<double>> DM, int n);	//函数4：评价函数

	vector<double> bestS(vector<double> Eval, int Length);	//函数5：搜索范围内最优评价值及其相应的位置函数

	vector<int> RandPosition(int n);		//函数6：产生Sharking操作位置函数
	vector<int> Swap(vector<int> S, vector<int> RP);	//函数7：交换算子
	vector<int> Flip(vector<int> S, vector<int> RP);	//函数8：翻转算子
	vector<int> Insert(vector<int> S, vector<int> RP);	//函数9：插入算子
};