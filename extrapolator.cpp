#include<iostream>
#include<vector>
#include"Gauss.h"

#pragma -O3

using namespace std;

long double summ_x (int first, int last, int back, vector<long double> x, vector<long double> mask) {
	long double result = 0;
	long double tmp = max(first, last);
	first = min(first, last);
	last = tmp;
	for (int i = first; i < last; i ++) {
		result += x[i - back] * mask[i];
	}
	return result;
}

long double summ_two_x (int first, int last, int back1, int back2, vector<long double> x, vector<long double> mask) {
	long double result = 0;
	long double tmp = max(first, last);
	first = min(first, last);
	last = tmp;
	for (int i = first; i < last; i ++) {
		result += x[i - back1] * x[i - back2] * mask[i];
	}
	return result;
}

vector<long double> write_and_solve_equation_for_f (vector<long double> x, int k, vector<long double> mask) {
	vector<vector<long double>> eqs = vector<vector<long double>> (k + 1);
	eqs[0].push_back(1);
	for (int i = 1; i <= k; i ++) {
		eqs[0].push_back(summ_x(x.size(), x.size() - k, i, x, mask));	
	}
	eqs[0].push_back(- summ_x(x.size(), x.size() - k, 0, x, mask));

	for (int i = 1; i < k + 1; i ++) {
		eqs[i].push_back(summ_x(x.size(), x.size() - k, i, x, mask));
		for (int j = 1; j <= k; j ++) {
			eqs[i].push_back(summ_two_x(x.size(), x.size() - k, i, j, x, mask));	
		}
		eqs[i].push_back(- summ_two_x(x.size(), x.size() - k, 0, i, x, mask));
	}

	return gaus_solve(eqs);
}

vector<long double> suggest(int distance, vector<long double> x, vector<long double> f) {
	vector<long double> suggestion = vector<long double> (distance, - f[0]);
	for (int i = 0; i < distance; i ++) {
		for (int j = 1; j < f.size(); j ++) {
			suggestion[i] -= f[j] * x[x.size() - j];
		}
		x.push_back(suggestion[i]);
	}
	return suggestion;
}

vector<long double> create_mask(int length, int important_length) {
	vector<long double> result = vector<long double>(length, 1);
	for (int i = important_length; i < length; i ++) {
		result[i] = important_length * 1.0 / i;
	}
	return vector<long double>(length, 1);
}

int main() {
//	vector<long double> x = {1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89};

//	vector<long double> x = {1, 1, 2, 3, 5, 8, 13, 21, 34, 55};
//	vector<long double> x = {59.0, 62.1, 65.2, 69.3, 70.9, 70.5, 70.7, 66.5, 65.8, 65.9, 66.4, 65.7, 64.1, 64.1, 62.1, 61.8, 62.3, 61.1, 61.0};
/*	vector<long double> x = {
72.7263,
73.6894,
74.5706,
74.4960,
74.7163,
75.1290,
77.0416,
76.2562,
74.6657,
73.9441};*/
	vector<long double> x = {3, 3, 6, 7, 10, 14, 17, 20, 20, 28, 34, 45, 59, 63, 93, 114, 147, 199, 253, 306, 367, 495, 658, 840, 1036, 1264, 1534, 1837, 2337, 2777, 3548, 4149, 4731, 5389, 6343, 7497, 8672, 10131, 11917, 13584, 15770, 18328, 21102, 27938, 32084, 36932, 42983, 47302, 52937, 57999, 62773, 68622, 74588, 80949, 87147, 93558, 99399, 106498, 114431, 124054, 134687, 145268, 155370, 165929, 177160, 187859, 198676, 209688, 221344
	//};
	, 232243}; 

//	cout << x.size() << endl;
	vector<long double> mask = create_mask(x.size(), 10);

	vector<long double> res = write_and_solve_equation_for_f (x, 5, mask);

	for (int i = 0; i < res.size(); i ++) {
		cout << res[i] << endl;
	}

	vector<long double> next = suggest(1, x, res);
for (int i = 0; i < next.size(); i ++) {
	cout << next[i] << endl;
}

	return 0;
}
