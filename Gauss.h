#ifndef GAUSS
#define GAUSS

#include<vector>

using namespace std;

vector<vector<long double>> exchange(vector<vector<long double>> eqs, int f, int s);

vector<vector<long double>> nullarize(vector<vector<long double>> eqs, int nullarizator, int nullarizated_point, int nullarizated_line);

vector<vector<long double>> triangalize(vector<vector<long double>> eqs);

long double summarize(vector<vector<long double>> eqs, int sumline, int begin_num, vector<long double> results);

vector<long double> gaus_solve(vector<vector<long double>> eqs);

#endif
