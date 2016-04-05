#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;


class solution{
public:
	double total_f;
	bool feasible;
	vector<int> trace;//the index of route which this request is assigned to


	vector<double> f;
	vector<vector<int>> route;

};

class raw{
public:
	int id;
	double x;
	double y;
	int d;
	int q;
	int e;
	int l;

	raw(double* a){
		id = (int)*a++;
		x = *a++;
		y = *a++;
		d = (int)*a++;
		q = (int)*a++;
		e = (int)*a++;
		l = (int)*a++;
	}
};

class tabu{
private:
	//raw input
	vector<raw> input;
	int m;
	int n;
	int T;
	int Q;
	int L;
	//processed data

	vector<vector<double>> t;//travel time
	vector<vector<double>> c;//distance

	//para
	double alpha;
	double beta;
	double gamma;
	double tau;
	double delta;
	
	double f;



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
			initial_dist();
	}
	void initial_dist(){
		for (int i = 0; i < 2 * n + 1;i++){
			vector<double> tem;
			for (int j = 0; j < 2 * n + 1;j++){
				double te = sqrt(pow((input[i].x - input[j].x), 2) + pow((input[i].y - input[j].y), 2));
				tem.push_back(te);
			}
			t.push_back(tem);
			c.push_back(tem);
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
			n = te[1]/2;
			T = te[2];
			Q = te[3];
			L = te[4];
		}
		else
			cout << "err in opening file";
		
	}
	//generate m routes ***  
	solution gen_init_ramdon(){
		solution s = solution();
		s.trace.push_back(-1);
		for (int i = 0; i < m ;i++){
			vector<int> tem;
			tem.push_back(0);
			s.route.push_back(tem);
			s.f.push_back(0);
		}
		srand((int)time(0));
		for (int i = 1; i <=n;i++)
		{
			int randnum = int(rand())%m;
			s.trace.push_back(randnum);
			s.route[randnum].push_back(i);
			s.route[randnum].push_back(i+n);
		}
		for (int i = 0; i < m;i++)
		{
			s.route[i].push_back(0);
		}
		return s;
	}

};

void main(){

	
	tabu a = tabu();
	solution t = a.gen_init_ramdon();
	cout << 1;
	
	}