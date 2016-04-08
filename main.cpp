#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;
template<class T>
int indexof(vector<T> a, T b){
	for (int i = 0; i < a.size(); i++){
		if (a[i]==b)
		{
			return i;
		}
	}
	return -1;
}
template<class T>
T max(T a, T b){
	if (a>b)
	{
		return a;
	}
	else return b;
}
template<class T>
T min(T a, T b){
	if (a > b)
	{
		return b;
	}
	else return a;
}
template<class T>
bool is_zero(T x){
	return x < 0.00001;
}

class route{
public:
	vector<int> T;
	vector<int> e;
	vector<int> l;
	vector<int> A;
	vector<int> B;
	vector<int> D;
	vector<int> W;
	bool fea;

	int duration;
	double f;
	int q;
	int d;
	int w;
	int t;
	route(){
		f = 0;
		q = 0;
		d = 0;
		w = 0;
		t = 0;
		duration = 0;
	}
	~route(){
		T.clear();
		e.clear();
		l.clear();
	}
};

class solution{
public:
	double total_f;
	vector<bool> fea;
	vector<int> trace;//the index of route which this request is assigned to

	vector<double> f;
	vector<vector<int>> rout;

	

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
			s.rout.push_back(tem);
			s.f.push_back(0);
		}
		srand((int)time(0));
		for (int i = 1; i <=n;i++)
		{
			int randnum = int(rand())%m;
			s.trace.push_back(randnum);
			s.rout[randnum].push_back(i);
			s.rout[randnum].push_back(i+n);
		}
		for (int i = 0; i < m;i++)
		{
			s.rout[i].push_back(0);
		}
		return s;
	}

	double objective_f(vector<int> a){
		route ret = route_evaluation(a);
		return ret.f + alpha*ret.q + beta*ret.d + gamma*ret.w + tau*ret.t;
	}

	double objective_f(route ret){
		return ret.f + alpha*ret.q + beta*ret.d + gamma*ret.w + tau*ret.t;
	}

	route route_evaluation(vector<int> a){
		route ret = route();
		for (int i = 0; i < a.size(); i++)
		{
			ret.e.push_back(input[a[i]].e);
			ret.l.push_back(input[a[i]].l);
			ret.T.push_back(0);
		}
		int QK = input[a[0]].q;
		
		vector<int> _D(2 * n + 1);
		vector<int> A(a.size());
		vector<int> B(a.size());
		vector<int> D(a.size());
		vector<int> W(a.size());

		A[0] = input[a[0]].e;
		B[0] = A[0];
		D[0] = B[0] + input[a[0]].d;
		W[0] = 0;
		
		for (int i = 1; i < a.size();i++)
		{
			A[i] = D[i - 1] + t[a[i - 1]][a[i]];
			B[i] = max(A[i], input[a[i]].e);
			D[i] = B[i] + input[a[i]].d;
			W[i] = input[a[i]].e - A[i];
			if (W[i]<0){
				W[i] = 0;
			}
			if (a[i] > 0 && a[i] <= n){
				_D[a[i]] = D[i];
			}
			if (a[i]>n){
				ret.T[i] = B[i] - _D[a[i] - n];
			}
		}
		//for 循环为debug内部
		for (int j = 0; j < a.size(); j++){
			int F0 = -1;
			int WSum = 0;
			for (int i = j; i < a.size(); i++){
				WSum += W[i];
				int _F0 = input[a[i]].l - B[i];
				if (a[i]>n){
					int _WSum = 0;
					for (int k = i; k >= j&&a[k] + n != a[i]; k--){
						_WSum += W[k];
						//}
						if (a[k] + n != a[i]){
							_WSum = 0;
							_F0 = min(_F0, L - ret.T[i] + _WSum);
						}
					}
				}
					if (i == a.size() - 1){
						_F0 = 0;
					}
					_F0 = WSum + max(0, _F0);
					if (F0<0 || F0>_F0){
						F0 = _F0;
					}
			}
			A[j] += F0;
			B[j] = max(A[j], input[a[j]].e);
			D[j] = B[j] + input[a[j]].d;
			W[j] = max(0, input[a[j]].e - A[j]);
			if (a[j] > 0 && a[j] <= n){
				_D[a[j]] = D[j];
			}
			for (int i = j + 1; i < a.size(); i++){
				A[i] = D[i - 1] + t[a[i - 1]][a[i]];
				B[i] = max(A[i], input[a[i]].e);
				D[i] = B[i] + input[a[i]].d;
				W[i] = input[a[i]].e - A[i];
				if (W[i] < 0){
					W[i] = 0;
				}
				if (a[i] > 0 && a[i] <= n){
					_D[a[i]] = D[i];
				}
				if (a[i] > n){
					ret.T[i] = B[i] - _D[a[i] - n];
				}
			}
		}
		//debug 结束
		for (int i = 1; i < a.size(); i++){
			ret.f += c[a[i - 1]][a[i]];
			
			QK += input[a[i]].q;
			if (QK > Q){
				ret.q += QK - Q;
			}
			if (B[i] > input[a[i]].l){
				ret.w += B[i] - input[a[i]].l;
			}
			if (a[i] > 0 && a[i] <= n){
				_D[a[i]] = D[i];
			}
			if (a[i] > n){
				ret.T[i] = B[i] - _D[a[i] - n];
				if (ret.T[i] > L){
					ret.t += ret.T[i] - L;
				}
			}
		}
		ret.duration = B[a.size() - 1] - B[0];
		if (ret.duration > T){
			ret.d = ret.duration - T;
		}
		ret.A = A;
		ret.B = B;
		ret.D = D;
		ret.W = W;


		return ret;
	}


	void simple_insertion_1(vector<int> rout, int v){
		int va, vb;
		int left = 1;
		int right = rout.size() - 1;

		va = indexof(rout, v - n);
		vb = indexof(rout, v + n);
		if (v > n && va > -1)
		{
			left = va + 1;
		}
		if (v <= n && vb > -1)
		{
			right = vb;
		}
		int m = 0, pos;
		for (int i = left; i <= right; i++)
		{
			vector<int> te(rout.size());
			copy(rout.begin(), rout.end(), te);
			te.insert(te.begin() + i, v);
			double f = objective_f(te);
			if (m > f){
				m = f;
				pos = i;
			}
		}
		rout.insert(rout.begin() + pos, v);
	}

	void simple_insertion_2(vector<int> rout, int v){
		int va, vb;
		if (input[v + n].l - input[v + n].e > input[v].l - input[v].e){
			va = v;
			vb = v + n;
		}
		else
		{
			vb = v;
			va = v + n;
		}
		simple_insertion_1(rout, va);
		simple_insertion_1(rout, vb);
	}

	bool is_fea(solution s){
		for (int i = 0; i < s.rout.size(); i++){
			if (!s.fea[i]){
				return false;
			}
		}
		return true;
	}

	bool is_fea_ret(route r){
		return is_zero(r.q) && is_zero(r.d) && is_zero(r.w) && is_zero(r.t);
	}

	double total_f(solution s){
		double f = 0;
		for (int i = 0; i < s.rout.size(); i++){
			route r = route_evaluation(s.rout[i]);
			s.fea.push_back(is_fea_ret(r));
			s.f.push_back(objective_f(r));
			f += s.f[i];
		}
		return f;
	}

	void solve(){
		alpha = 10;
		beta = 0.1;
		gamma = 1;
		tau = 1;
		delta = 0.1;
		f = 0;
		solution s = gen_init_ramdon();

		int iter_best, iter_end = 400;
		int iter_internal = 10;
		
	}

};

void main(){

	
	tabu a = tabu();
	solution t = a.gen_init_ramdon();
	cout << 1;
	/*
	vector<int> myints = { 10, 20, 30, 40 };
	int i;
	i = indexof(myints, 40);
	cout << i;
	system("pause");
	*/
}