#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;


struct solution{
	vector<double> f;
	vector<int> route;
};

class raw{
	int id;
	double x;
	double y;
	int d;
	int q;
	int e;
	int l;
public:
	raw(double* a){
		id = *a++;
		x = *a++;
		y = *a++;
		d = *a++;
		q = *a++;
		e = *a++;
		l = *a++;
	}
};

class tabu{
private:

	vector<raw> input;
	int m;
	int n;
	int T;
	int Q;
	int L;
public:
	tabu(void){
		ifstream in("pr01", ios::in);
		string inp;
		initial_tabu(in);
			for (int j = 0; j <=n;j++)
			{
				while (getline(in, inp))
				{
					stringstream stringin(inp);
					double *a = new double[7];
					for (int j = 0; j < 7;j++)
					{
						stringin >> a[j];
					}
					input.push_back(raw(a));
				}
			}
	}

	void initial_tabu(ifstream & a){
		if (a)
		{
			string input;
			getline(a, input);
			stringstream stringin(input);
			int tem, i = 0;
			int te[5];
			while (stringin >> tem)
			{
				te[i] = tem;
				i++;
			}
			m = te[0];
			n = te[1];
			T = te[2];
			Q = te[3];
			L = te[4];
		}
		else
			cout << "err in opening file";
		
	}
	void 

};

void main(){
	tabu a = tabu();
	cout << 1;
}