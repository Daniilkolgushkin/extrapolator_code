#include<vector>

#pragma -O3

using namespace std;

vector<vector<long double>> exchange(vector<vector<long double>> eqs, int f, int s) {
	vector<long double> tmp = eqs[f];
	eqs[f] = eqs[s];
	eqs[s] = tmp;
	return eqs;
}

vector<vector<long double>> nullarize(vector<vector<long double>> eqs, int nullarizator, int nullarizated_point, int nullarizated_line) {
	long double nullarization_koefficient = eqs[nullarizated_line][nullarizated_point] / eqs[nullarizator][nullarizated_point];
	for (int i = nullarizated_point; i < eqs[nullarizated_line].size(); i ++) {
		eqs[nullarizated_line][i] = eqs[nullarizated_line][i] - nullarization_koefficient * eqs[nullarizator][i];
	}
	return eqs;
}

vector<vector<long double>> triangalize(vector<vector<long double>> eqs) {
	for (int i = 0; i < eqs.size(); i ++) {
		if (eqs[i][i] == 0) {
			for (int j = i + 1; j < eqs.size(); j ++) {
				if (eqs[j][i] != 0) {
					eqs = exchange(eqs, i, j);
				}
			}
		}
		for (int j = i + 1; j < eqs.size(); j ++) {
			eqs = nullarize(eqs, i, i, j);
		}
	}
	return eqs;
}

long double summarize(vector<vector<long double>> eqs, int sumline, int begin_num, vector<long double> results) {
	long double result = 0;
	for (int i = begin_num; i < eqs[sumline].size() - 1; i ++) {
		result += eqs[sumline][i] * results[i];
	}
	return result;
}

vector<long double> gaus_solve(vector<vector<long double>> eqs) {
	eqs = triangalize(eqs);
	vector<long double> result = vector<long double>(eqs[0].size() - 1);
	for (int i = eqs.size() - 1; i >= 0; i --) {
		result[i] = (eqs[i][eqs[i].size() - 1] - summarize(eqs, i, i + 1, result)) / (eqs[i][i]);
	}
	return result;
}

/*vector<vector<long double>> rectangalize(vector<vector<long double>> eqs) {
	int maximal_length = 0;
	for (int i = 0; i < eqs.size(); i ++) {
		if (eqs[i].size() > maximal_length) {
			maximal_length = eqs[i].size();
		}
	}
	for (int i = 0; i < eqs.size(); i ++) {
		for (int j = eqs[i].size(); j < maximal_length; j ++) {
			eqs[i].push_back(0);
		}
	}
	return eqs;
}

vector<vector<long double>> read(string input_filename) {
	ifstream input(input_filename);
	string current_number;
	vector<vector<long double>> eqs = vector<vector<long double>>(1);
	while(!input.eof()) {
		input >> current_number;
		if (current_number == ";") {
			eqs.push_back({});

		}
		else if (current_number[current_number.size() - 1] == ';'){
			current_number = current_number.substr(0, current_number.size() - 1);
			eqs[eqs.size() - 1].push_back(stod(current_number));
			eqs.push_back({});
		}
		else {
			eqs[eqs.size() - 1].push_back(stod(current_number));

		}
	}
	eqs[eqs.size() - 1].pop_back();
	return rectangalize(eqs);
}

void write(string output_filename, vector<long double> result) {
	ofstream output(output_filename);
	for (int i = 0; i < result.size(); i ++){
		output << result[i] << endl;
	}
}

int main() {
	string input_filename;
	string output_filename;
	cin >> input_filename;
	cin >> output_filename;
	vector<vector<long double>> in = read(input_filename);
	vector<long double> res = gaus_solve(in);
	write(output_filename, res);
	return 0;
}*/
