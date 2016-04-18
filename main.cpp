#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;
clock_t start_t, end_t;


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

double dmax(double a, int b){
	if (a > b)
	{
		return a;
	}
	else return (double)b;
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

double sum_vec(vector<double> a, int x){
	double res = 0;
	for (vector<double>::iterator it = a.begin()+x-1; it!= a.end()-1; it++){
		res += *it;
	}
	return res;
}


class route{
public:
	vector<double> T;
	vector<double> e;
	vector<double> l;
	vector<double> A;
	vector<double> B;
	vector<double> D;
	vector<double> W;
	vector<double> L;
	vector<double> P;

	bool fea;

	int duration;
	double f;
	double q;
	double d;
	double w;
	double t;
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

	string nam;


	tabu(string fi){
		nam = fi;
		ifstream in(fi, ios::in);
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
		route ret;
		int size = a.size();
		/*for (size_t i = 0; i < size; i++)
		{
			ret.e.push_back(input[a[i]].e);
			ret.l.push_back(input[a[i]].l);
			ret.T.push_back(0);
		}*/
		ret.T.resize(size);
		int QK = input[a[0]].q;
		
		vector<int> _D(2 * n + 1);
		vector<double> A(size);
		vector<double> B(size);
		vector<double> D(size);
		vector<double> W(size);
		

		A[0] = input[a[0]].e;
		B[0] = A[0];
		D[0] = B[0] + input[a[0]].d;
		W[0] = 0;
		int count = 1;
		vector<int>::iterator ia = a.begin() + 1;

		for (vector<double>::iterator ita = A.begin() + 1, itb = B.begin() + 1, itd = D.begin() + 1, itw = W.begin() + 1; ita != A.end(); ita++, itb++, itd++, itw++)
		{

			*ita = *(itd-1) + t[*(ia-1)][*ia];
			
			*itb = dmax(*ita, input[*ia].e);
			*itd = *itb + input[*ia].d;
			*itw = input[*ia].e - *ita;
			if (*itw<0){
				*itw = 0;
			}
			if (*ia > 0 && *ia <= n){
				_D[*ia] = *itd;
			}
			if (*ia>n){
				ret.T[count] = *itb - _D[*ia - n];
			}
			count++;
			ia++;
		}
		ia = a.begin();
		count = 1;
		//for 循环为debug内部

		vector<double>::iterator itw = W.begin(), ita = A.begin(), itb = B.begin(), itd = D.begin();


		for (int j = 0; j < size; j++){
			double F0 = -1;
			int WSum = 0;
			for (int i = j; i < size; i++){
				WSum += W[i];
				double _F0 = input[a[i]].l - B[i];
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
					if (i == size - 1){
						_F0 = 0;
					}
					_F0 = WSum + dmax(0, _F0);
					if (F0<0 || F0>_F0){
						F0 = _F0;
					}
			}

			*ita += F0;
			//*itb = *itb;
			*itb = dmax(*ita, input[*ia].e);
			//*itb += F0;
			*itd = *itb + input[*ia].d;
			*itw = dmax(0, input[*ia].e - *ita);
			if (*ia > 0 && *ia <= n){
				_D[*ia] = *itd;
			}

			
			vector<double>::iterator itw1 = W.begin()+j+1, ita1 = A.begin()+j+1, itb1 = B.begin()+j+1, itd1 = D.begin()+j+1;
			vector<int>::iterator iatp1 = a.begin() + j + 1;
			for (int i = j + 1; i < size; i++){
				*ita1 = *(itd1-1) + t[*(iatp1-1)][*iatp1];
				*itb1 = dmax(*ita1, input[*iatp1].e);
				*itd1 = *itb1 + input[*iatp1].d;
				*itw1 = input[*iatp1].e - *ita1;
				if (*itw1 < 0){
					*itw1 = 0;
				}
				if (*iatp1 > 0 && *iatp1 <= n){
					_D[*iatp1] = *itd1;
				}
				if (*iatp1 > n){
					ret.T[i] = *itb1 - _D[*iatp1 - n];
				}

				ita1++;
				itb1++;
				itd1++;
				itw1++;
				iatp1++;
			}
			
			/*
			
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
			*/
			


			ita++;
			itb++;
			itd++;
			itw++;
			ia++;



		}
		B[0] = B[1] - t[a[0]][a[1]];
		A[0] = B[0];
		D[0] = B[0];
		//debug 结束
		ia = a.begin() + 1;
		for (vector<double>::iterator itb = B.begin()+1,itt = ret.T.begin()+1,itd=D.begin()+1; itb!=B.end(); itb++,itt++,itd++){
			
			ret.f += c[*(ia-1)][*ia];
			
			QK += input[*ia].q;
			if (QK > Q){
				ret.q += QK - Q;
			}
			if (*itb > input[*ia].l){
				ret.w += *itb - input[*ia].l;
			}
			if (*ia > 0 && *ia <= n){
				_D[*ia] = *itd;
			}
			if (*ia > n){
				*itt = *itb - _D[*ia - n];
				if (*itt > L){
					ret.t += *itt - L;
				}
			}
			ia++;
		}
		ret.duration = B[size - 1] - B[0];
		if (ret.duration > T){
			ret.d = ret.duration - T;
		}
		ret.A = A;
		ret.B = B;
		ret.D = D;
		ret.W = W;


		return ret;
	}

	route route_evaluation_new(vector<int> a){
		route ret;

		for (size_t i = 0; i < a.size(); i++)
		{
			ret.e.push_back(input[a[i]].e);
			ret.l.push_back(input[a[i]].l);
		}
		vector<double> _D(2 * n + 1);

		ret.A.resize(a.size());
		ret.B.resize(a.size());
		ret.D.resize(a.size());
		ret.W.resize(a.size());
		ret.P.resize(a.size());

		ret.A[0] = input[a[0]].e;
		ret.B[0] = ret.A[0];
		ret.D[0] = ret.B[0] + input[a[0]].d;
		ret.W[0] = 0;
		ret.P[0] = 0;
		for (size_t i = 1; i < a.size(); i++){
			ret.A[i] = ret.D[i - 1] + t[a[i - 1]][a[i]];
			ret.B[i] = dmax(ret.A[i], input[a[i]].e);
			ret.D[i] = ret.B[i] + input[a[i]].d;
			_D[a[i]] = ret.D[i];
			ret.W[i] = ret.B[i] - ret.A[i];
		}
		for (size_t i = 1; i < a.size(); i++){
			if (a[i] > n || a[i] == 0){
				ret.P[i] = 0;
			}
			else{
				ret.P[i] = -ret.B[i] + _D[a[i] + n];
			}
		}

		ret.D[0] = input[a[0]].e + dmin(F_new(0, a.size()-1, a, ret), sum_vec(ret.W, 1));
		for (size_t i = 1; i < a.size(); i++){
			ret.A[i] = ret.D[i - 1] + t[a[i - 1]][a[i]];
			ret.B[i] = dmax(ret.A[i], input[a[i]].e);
			ret.D[i] = ret.B[i] + input[a[i]].d;
			_D[a[i]] = ret.D[i];
			ret.W[i] = ret.B[i] - ret.A[i];
		}

		for (size_t i = 1; i < a.size(); i++){
			if (0 < a[i] && a[i] <= n){
				ret.P[i] = -ret.B[i] + _D[a[i] + n];
			}
		}

		for (size_t i = 1; i < a.size()-1; i++){
			if (a[i]>n||a[i]==0){ continue; }
			double Fi = F_new(i, a.size()-1, a, ret);
			ret.B[i] += dmin(Fi, sum_vec(ret.W, i + 1));
			ret.D[i] = ret.B[i] + input[a[i]].d;
			_D[a[i]] = ret.D[i];
			ret.P[i] = -ret.B[i] + _D[a[i] + n];
			for (int j = i + 1; j < a.size(); j++){
				ret.A[j] = ret.D[j - 1] + t[a[j - 1]][a[j]];
				ret.B[j] = dmax(ret.A[j], input[a[j]].e);
				ret.D[j] = ret.B[j] + input[a[j]].d;
				_D[a[j]] = ret.D[j];
				ret.W[j] = ret.B[j] - ret.A[j];
			}
			for (int j = i + 1; j < a.size(); j++){
				if (0 < a[i] && a[i] <= n){
					ret.P[i] = -ret.B[i] + _D[a[i] + n];
				}
			}

		}
		double vq = 0, vd = 0, vw = 0, vt = 0;
		vector<double> temp(a.size());
		for (int i = 0; i < a.size(); i++){
			double _Q = 0;
			_Q += input[a[i]].q;
			if (_Q > Q){
				vq += _Q - Q;
			}
			if (ret.B[i] > input[a[i]].l){
				vw += ret.B[i] - input[a[i]].l;
			}
			else if (ret.B[i] < input[a[i]].e){
				vw += -ret.B[i] + input[a[i]].e;
			}
			if (a[i] > n){
				double ti = ret.B[i] - _D[a[i] - n];
				if (ti > L){
					vt += ti - L;
				}
			}
		}
		ret.duration = ret.D[a.size() - 1] - ret.D[0];
		if (ret.duration > T){
			vd += ret.duration - T;
		}
		ret.q = vq;
		ret.d = vd;
		ret.w = vw;
		ret.t = vt;


		return ret;
	}

	double F(int a, int q, vector<int> r,route ret){
		double res = DBL_MAX;
		for (int i = a; i <= q; i++){
			double tem2 = 0;
			for (int j = a; j < i; j++){
				tem2 += t[r[j]][r[j + 1]];
			}
			double tem = input[r[i]].l - (ret.B[a] + tem2);
			if (res > tem){
				res = tem;
			}
		}
		return res;
	}

	double F_new(int a, int q, vector<int> r, route ret){
		double res = DBL_MAX;
		for (int i = a; i <= q; i++){
			double tem2 = 0;
			for (int j = a+1; j <= i; j++){
				tem2 += ret.W[j];
			}
			double tem = tem2 + dmax(dmin(input[r[i]].l-ret.B[i],L-ret.P[i]),0);
			if (res > tem){
				res = tem;
			}
		}
		return res;
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
	double c_f(solution &s){
		double f = 0;
		int si = s.rout.size();
		for (int i = 0; i < si; i++){
			route r = route_evaluation(s.rout[i]);
			f += r.f;
		}
		return f;
	}

	double tsolve(){

		//ofstream out;
		//out.open(nam + "_solve",ios::app);
		//vector<double> res(10);
		alpha = 1;
		beta = 0.1;
		gamma = 10;
		tau = 1;
		delta = 0.5;
		f = 0;
		solution s = gen_init_ramdon();

		start_t = clock();

		int iter_best, iter_end = 500;
		int iter_internal = 10;
		int theta = (int)7.5*log10(iter_end);
		int lambda = 0.015;
		for (int i = 0; i < n+1; i++){
			vector<int> tem(m);
			for (int j = 0; j < m; j++){
				tem[j] = 0;
			}
			rho.push_back(tem);
		}

		for (int iter = 0; iter < iter_end; iter++){
			//delta = 7.5*log10(iter);


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
					newf = f01 + f11 - f00 - f10 + rho[i][j] * lambda*total_f(s)*sqrt(n*m);
					//cout << "i:" << i << " j:" << j << " newf: " << newf << " minf:" << minf << endl;
					//cout << " f01:" << f01 << " f11:" << f11 << " f00" << f00 << " f10" << f10 << " rho:" << rho[i][j] << endl;
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
				//cout << "internal iter"<<endl;
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
			
			//out << "times: " << iter << " feasible: " << fea << " f: " << f << endl;
			
			if (fea){
				ratio = 1 / ratio;
				if (!bestfound || bestf>f){
					bestfound = true;
					bestf = f;
					bests = s;
					//time_best = time();

					iter_best = iter;
					//cout << "better f" << f << endl;

				}
			}
			alpha *= ratio;
			beta *= ratio;
			gamma *= ratio;
			tau *= ratio;
		}
		end_t = clock();
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
			return bestf;
			//out << "f:" << bestf << endl;
			//out << "duration:" << dura << endl;
			//out << "waiting:" << waiting << endl;
			//out << "transit:" << transit << endl;
		}
		//out.close();
	}

	double dsum(vector<double> a){
		double res = 0;
		for (vector<double>::iterator it = a.begin(); it != a.end(); it++){
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
	//tabu a("pr02");
	//route t = a.route_evaluation(vector < int > {0,10,11,35,34,0});
	
	double res = 0;
	for (int j = 19; j < 20; j++){
		
		for (int i = 1; i < 2; i++){
			
			string s = "pr";
			if (j < 10){ s = s + "0" + to_string(j); }
			else{ s += to_string(j); }
			tabu a(s);

			double r = a.tsolve();
			res += r;
			double time = (double)(end_t - start_t) / CLOCKS_PER_SEC / 60;
			//ofstream out;
			//out.open(a.nam + "_solve", ios::app);
			cout << "time(min): " << time <<" f: "<<r<< endl;
			//out.close();
		}
	}
	//cout << res  << endl;
	
	/*
	route r = a.route_evaluation(vector < int > {0, 10, 11, 35, 34, 0});
	double res = 0;
	for (vector<double>::iterator it = r.W.begin(); it != r.W.end(); it++){
		cout << *it << endl;
		res += *it;
	}
	cout << "total W" << res << endl;
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