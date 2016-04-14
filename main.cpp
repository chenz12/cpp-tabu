#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;


class dpair{
public:
	int x;
	int y;
	dpair(int a, int b){
		x = a;
		y = b;
	}

};
template<class T>
int tindexof(vector<T> a, T b){
	for (int i = 0; i < a.size(); i++){
		if (equal_pair(a[i],b))
		{
			return i;
		}
	}
	return -1;
}

bool equal_pair(dpair a, dpair b){
	if (a.x == b.x&&a.y == b.y){
		return true;
	}
	return false;
}
bool equal_pair(int a, int b){
	if (a==b){
		return true;
	}
	return false;
}

template<class T>
T dmax(T a, T b){
	if (a>b)
	{
		return a;
	}
	else return b;
}
template<class T>
T dmin(T a, T b){
	if (a > b)
	{
		return b;
	}
	else return a;
}
template<class T>
bool tis_zero(T x){
	return x < 0.00001;
}

vector<int> tremove(vector<int> &a, int i) {
	for (vector<int>::iterator it = a.begin(); it != a.end(); it++) {
		if (*it == i) {
			it = a.erase(it);
		}
	}
	return a;
}

template<class T>
vector<T> routecopy(vector<T> a) {
	vector<T> res;
	for (vector<T>::iterator it = a.begin(); it != a.end(); it++) {
		res.push_back(*it);
	}
	return res;
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
public:
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

	//tabu list
	vector<dpair> tabulist;

	//rho
	vector<vector<int>> rho;

	//process parameter
	bool bestfound = false;
	double bestf;
	solution bests;




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
			s.fea.push_back(false);
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
		for (size_t i = 0; i < a.size(); i++)
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
		
		for (size_t i = 1; i < a.size();i++)
		{
			A[i] = D[i - 1] + t[a[i - 1]][a[i]];
			B[i] = dmax(A[i], input[a[i]].e);
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
							_F0 = dmin(_F0, L - ret.T[i] + _WSum);
						}
					}
				}
					if (i == a.size() - 1){
						_F0 = 0;
					}
					_F0 = WSum + dmax(0, _F0);
					if (F0<0 || F0>_F0){
						F0 = _F0;
					}
			}
			A[j] += F0;
			B[j] = dmax(A[j], input[a[j]].e);
			D[j] = B[j] + input[a[j]].d;
			W[j] = dmax(0, input[a[j]].e - A[j]);
			if (a[j] > 0 && a[j] <= n){
				_D[a[j]] = D[j];
			}
			for (int i = j + 1; i < a.size(); i++){
				A[i] = D[i - 1] + t[a[i - 1]][a[i]];
				B[i] = dmax(A[i], input[a[i]].e);
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


	void simple_insertion_1(vector<int> &rout, int v){
		int va, vb;
		int left = 1;
		int right = rout.size() - 1;

		va = tindexof(rout, v - n);
		vb = tindexof(rout, v + n);
		if (v > n && va > -1)
		{
			left = va + 1;
		}
		if (v <= n && vb > -1)
		{
			right = vb;
		}
		int m = 0, pos;
		double f;
		for (int i = left; i <= right; i++)
		{
			vector<int> te = rout;
			te.insert(te.begin() + i, v);
			f = objective_f(te);
			
			if (m > f||m==0){
				m = f;
				pos = i;
			}
		}
		rout.insert(rout.begin() + pos, v);
	}

	void simple_insertion_2(vector<int> &rout, int v){
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
		return tis_zero(r.q) && tis_zero(r.d) && tis_zero(r.w) && tis_zero(r.t);
	}

	double total_f(solution &s){
		double f = 0;
		for (int i = 0; i < s.rout.size(); i++){
			route r = route_evaluation(s.rout[i]);
			s.fea[i]=is_fea_ret(r);
			s.f[i]=objective_f(r);
			f += s.f[i];
		}
		s.total_f = f;
		return f;
	}

	void tsolve(){
		alpha = 10;
		beta = 0.1;
		gamma = 1;
		tau = 1;
		delta = 0.1;
		f = 0;
		solution s = gen_init_ramdon();

		int iter_best, iter_end = 40;
		int iter_internal = 10;
		int theta = iter_end / 40;
		int lambda = 100;
		for (int i = 0; i < n+1; i++){
			vector<int> tem(m);
			for (int j = 0; j < m; j++){
				tem[j] = 0;
			}
			rho.push_back(tem);
		}

		for (int iter = 0; iter < iter_end; iter++){
			f = total_f(s);
			double minf = 10000;
			int mini = -1;
			int minj = -1;
			double newf;

			//bool dbg;

			for (int i = 1; i <= n; i++){
				for (int j = 0; j < m; j++){
					if (j == s.trace[i]){
						continue;
					}
					dpair ta = dpair(i, j);
					int k0 = s.trace[i];
					int k1 = j;
					vector<int>* s0 = &s.rout[k0];
					vector<int>* s1 = &s.rout[k1];
					double f00 = s.f[k0];
					double f10 = s.f[k1];
					vector<int> r0 = *s0;
					vector<int> r1 = *s1;
					tremove(r0, i);
					tremove(r0, i + n);
					simple_insertion_2(r1, i);
					double f01 = objective_f(r0);
					double f11 = objective_f(r1);
					newf = f01 + f11 - f00 - f10 + rho[i][j] * lambda;
					cout << "i:" << i << " j:" << j << " newf: " << newf << " minf:" << minf << endl;
					cout << " f01:" << f01 << " f11:" << f11 << " f00" << f00 << " f10" << f10 << " rho:" << rho[i][j] << endl;
					if (tindexof(tabulist, ta) > -1 && f + newf > bestf){
						//dbg = true;
						continue;
					}

					if (minf > newf){
						minf = newf;
						mini = i;
						minj = j;
					}
				}
			}
			if (mini == -1&&minj==-1){
				cout << "buging" << endl;
			//	continue;
			}
			int i = mini;
			int j = minj;
			rho[i][j]++;
			int k0 = s.trace[i];
			int k1 = j;
			vector<int>* s0 = &s.rout[k0];
			vector<int>* s1 = &s.rout[k1];
			tremove(*s0, i);
			tremove(*s0, i + n);
			simple_insertion_2(*s1, i);
			s.trace[i] = j;
			dpair ta = dpair(i, k0);
			if (tindexof(tabulist, ta) == -1){
				tabu_push(tabulist, ta, theta);
			}
			if (iter%iter_internal == 0){
				cout << "internal iter"<<endl;
				for (int k = 1; k <=n; k++){
					vector<int> *ss = &s.rout[s.trace[k]];
					tremove(*ss, k);
					tremove(*ss, k + n);
					simple_insertion_2(*ss, k);
				}
			}
			f = total_f(s);
			bool fea = is_fea(s);
			double ratio = 1 + delta;
			if (fea){
				ratio = 1 / ratio;
				if (!bestfound || bestf>f){
					bestfound = true;
					bestf = f;
					bests = s;
					//time_best = time();

					iter_best = iter;
					cout << "better f" << f << endl;

				}
			}
			alpha *= ratio;
			beta *= ratio;
			gamma *= ratio;
			tau *= ratio;
		}

		//time_end = time();
		if (bestfound){
			double dura = 0;
			double waiting = 0;
			double transit = 0;
			for (int i = 0; i < s.rout.size(); i++){
				route ret = route_evaluation(s.rout[i]);
				dura += ret.duration;
				waiting += dsum(ret.W);
				transit += dsum(ret.T);
			}
			cout << "f:" << bestf << endl;
			cout << "duration:" << dura << endl;
			cout << "waiting:" << waiting << endl;
			cout << "transit:" << transit << endl;
		}
	}

	int dsum(vector<int> a){
		int res = 0;
		for (vector<int>::iterator it = a.begin(); it != a.end(); it++){
			res += *it;
		}
		return res;
	}

	void tabu_push(vector<dpair> &a, dpair b, int c){
		if (a.size() < c){
			a.push_back(b);
		}
		else if (a.size() == c){
			a.erase(a.begin());
			a.push_back(b);
		}
		else {
			cout << "error";
		}
	}
};

void main(){

	
	tabu a = tabu();
	solution t = a.gen_init_ramdon();
	a.tsolve();
	/*
	for (int i = 0; i < 20; i++){
		cout << "cost1：" << a.objective_f(t.rout[0]) << endl;
		cout << "cost2：" << a.objective_f(a.route_evaluation(t.rout[0])) << endl;
		cout << "feasible: " << a.is_fea_ret(a.route_evaluation(t.rout[0])) << endl;
	}
	cout << 1;*/
	system("pause");
	/*
	vector<int> myints = { 10, 20, 30, 40 };
	int i;
	i = indexof(myints, 40);
	cout << i;
	system("pause");
	*/
}