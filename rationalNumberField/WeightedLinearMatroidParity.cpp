#include <iostream>
#include <cstdio>
#include <math.h>
#include <string>
#include <algorithm>
#include <functional>
#include <vector>
#include <set>
#include <fstream>
#include <map>
#include <queue>
#include <time.h>
#include <stack>
#define INF 9999999999
#define EPS 1.0e-10
#include "Field.h"
#include "RationalNumber.h"
#include "Matrix.h"
#include "Tree.h"

#include <boost/multiprecision/cpp_int.hpp>

using namespace std;

typedef long long ll;
typedef std::pair<ll, ll> pr;

ll blossom_number;
ll vertex_number;
ll ordering_number;

ll Blossom_num = 0;
ll Graft_num = 0;
ll Augment_num = 0;
ll Search_num = 0;
ll DualUpdate_num = 0;

num_type maximum_absolute_num = 0;
bool overflow = false;

//pivoting around a pair p
void Pivoting_around_p(Matrix &C, pr p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase)
{
	ll x = p.first;
	ll y = p.second; //x in B and y not in B
	//cout << "x:" << x << " y:" << y << "number:";
	//C.X[x][y].output();
	//cout << endl;

	/*
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << "Ba.size is " << Ba.size() << " NBa.size() is " << NBa.size() << endl;
	cout << Ba[x] << " " << NBa[y] << endl;

	cout << "C:" << endl;
	C.output_matrix();
	cout << endl;
	*/

	for (ll i = 0; i < C.row; i++)
	{
		for (ll j = 0; j < C.col; j++)
		{
			//cout << i << " " << j << endl;
			if (i != x && j != y)
			{
				C.X[i][j] = C.X[i][j] - C.X[i][y] * C.X[x][j] / C.X[x][y];
			}
		}
	}
	for (ll i = 0; i < C.row; i++)
	{
		if (i != x)
		{
			C.X[i][y] = -C.X[i][y] / C.X[x][y];
		}
	}
	for (ll j = 0; j < C.col; j++)
	{
		if (j != y)
		{
			C.X[x][j] = C.X[x][j] / C.X[x][y];
		}
	}
	C.X[x][y] = C.X[x][y].inverse();

	ll x_sub;
	ll y_sub;
	for (ll i = 0; i < Ba_sub.size(); i++)
	{
		if (Ba_sub[i] == Ba[x])
		{
			x_sub = i;
			break;
		}
	}
	for (ll i = 0; i < NBa_sub.size(); i++)
	{
		if (NBa_sub[i] == NBa[y])
		{
			y_sub = i;
			break;
		}
	}
	swap(Ba[x], NBa[y]);
	swap(Ba_sub[x_sub], NBa_sub[y_sub]);
	swap(inBase[Ba[x]], inBase[NBa[y]]);

	num_type C_mx = C.maximum_number();
	//cout << "the maximum absolute number in elements of C:";
	//cout << C_mx << endl;

	maximum_absolute_num = max(maximum_absolute_num, C_mx);
	/*
	cout << "C:" << endl;
	C.output_matrix();
	cout << endl;
	*/
}

//pivoting around a matching M
void Pivoting_around_M(Matrix &C, vector<pr> &M, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase)
{
	//cout << "-----pivoting around M start-----" << endl;

	////
	/*
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}

	cout << "C:" << endl;
	C.output_matrix();
	cout << endl;
	*/

	vector<bool> flag_Ba;
	vector<bool> flag_NBa;
	flag_Ba.resize(C.row, false);
	flag_NBa.resize(C.col, false);
	vector<ll> A;
	vector<ll> B;
	//C.output_matrix();
	//cout << endl;
	for (ll i = 0; i < M.size(); i++)
	{
		//cout << M[i].first << " " << M[i].second << endl;
		A.push_back(M[i].first);
		B.push_back(M[i].second);
		flag_Ba[M[i].first] = true;
		flag_NBa[M[i].second] = true;
	}
	vector<ll> NA;
	vector<ll> NB;
	for (ll i = 0; i < Ba.size(); i++)
	{
		if (!flag_Ba[i])
		{
			NA.push_back(i);
		}
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		if (!flag_NBa[i])
		{
			NB.push_back(i);
		}
	}

	Matrix alpha = C.sub_matrix(A, B);
	Matrix beta = C.sub_matrix(A, NB);
	Matrix gamma = C.sub_matrix(NA, B);
	Matrix delta = C.sub_matrix(NA, NB);

	/////
	/*
	cout << "alpha:" << endl;
	alpha.output_matrix();
	matrix_inverse(alpha).output_matrix();
	*/
	/////

	Matrix alpha_inverse = matrix_inverse(alpha);

	/////
	/*
	bool flag = 1;
	for(ll i=0;i<alpha_inverse.row;i++){
		for(ll j=0;j<alpha_inverse.col;j++){
			if(!alpha_inverse.X[i][j].isZero())
			{
				flag=0;break;
			}
		}
	}
	if(flag){
		return false;
	}
	*/
	/////

	alpha = alpha_inverse;
	Matrix beta_ = alpha_inverse * beta;
	Matrix gamma_ = -gamma * alpha_inverse;
	delta = delta - gamma * alpha_inverse * beta;

	for (ll i = 0; i < A.size(); i++)
	{
		for (ll j = 0; j < B.size(); j++)
		{
			C.X[A[i]][B[j]] = alpha.X[i][j];
		}
	}
	for (ll i = 0; i < A.size(); i++)
	{
		for (ll j = 0; j < NB.size(); j++)
		{
			C.X[A[i]][NB[j]] = beta_.X[i][j];
		}
	}
	for (ll i = 0; i < NA.size(); i++)
	{
		for (ll j = 0; j < B.size(); j++)
		{
			C.X[NA[i]][B[j]] = gamma_.X[i][j];
		}
	}
	for (ll i = 0; i < NA.size(); i++)
	{
		for (ll j = 0; j < NB.size(); j++)
		{
			C.X[NA[i]][NB[j]] = delta.X[i][j];
		}
	}
	for (ll i = 0; i < M.size(); i++)
	{
		swap(Ba[M[i].first], NBa[M[i].second]);
		swap(inBase[Ba[M[i].first]], inBase[NBa[M[i].second]]);
		//swap(Ba[A[i]],NBa[B[i]])
	}
	/*
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;
	cout << "C:" << endl;
	C.output_matrix();
	cout << endl;
	*/
	num_type C_mx = C.maximum_number();
	//cout << "the maximum absolute number in elements of C:";
	//cout << C_mx << endl;

	maximum_absolute_num = max(maximum_absolute_num, C_mx);

	//return true;
	//cout << "-----pivoting around M end------" << endl;
}

//judge if {i, j} is a line or not
bool isLine(ll i, ll j)
{
	if (abs(i - j) == 1 && min(i, j) % 2 == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//get the mate of i
ll mate(ll i)
{
	if (i % 2 == 0)
	{
		return i + 1;
	}
	else
	{
		return i - 1;
	}
}

//judge if the vertex i is a single vertex or not
bool isSingle(ll i, vector<ll> &K)
{
	if (K[i] == -i)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//judge if K(u) and K(v) matches each other
bool judge_K(ll u, ll v, vector<ll> &K)
{
	if (K[u] > 0 && K[v] > 0)
	{
		return K[u] == K[v];
	}
	else if (K[u] <= 0 && K[v] <= 0)
	{
		return isLine(-K[u], -K[v]);
	}
	else
	{
		return false;
	}
}

//judge if there exists a source line in a given vector
//Note: N is the number of columns of A.
bool SourceLine_vector(ll N, vector<ll> &Ba, vector<bool> &inBase)
{
	bool ret = false;
	for (ll i = 0; i < Ba.size(); i++)
	{
		if (Ba[i] < N && !inBase[mate(Ba[i])])
		{
			ret = true;
		}
	}
	return ret;
}

//judge if there exists a source line in a given set
bool SourceLine_set(ll N, set<ll> &St, vector<bool> &inBase)
{
	bool ret = false;
	for (auto itr = St.begin(); itr != St.end(); itr++)
	{
		if (*itr < N && inBase[mate(*itr)] != inBase[*itr])
		{
			ret = true;
		}
	}
	return ret;
}

//get the path from v to a source vertex
vector<ll> Path(ll v, vector<vector<ll>> &P)
{
	//cout << "-----Path start-----" << endl;
	/* cout << "v is " << v << endl;
	cout << "P:" << endl;
	for (ll i = 0; i < P.size(); i++)
	{
		cout << "P[" << i << "]:";
		for (ll j = 0; j < P[i].size(); j++)
		{
			cout << P[i][j] << " ";
		}
		cout << endl;
	} */

	vector<ll> path;
	path.push_back(v);
	ll cur = v;
	while (P[cur][0] != cur)
	{
		for (ll i = P[cur].size() - 1; i >= 0; i--)
		{
			path.push_back(P[cur][i]);
		}
		cur = P[cur][0];
	}

	//cout << "-----path end-----" << endl;

	return path;
}

bool comp(const pr &a, const pr &b)
{
	return a.first > b.first;
}

//add new vertex to V
ll new_vertex(vector<ll> &K, vector<vector<ll>> &P, Matrix &C, Matrix &Q, vector<ll> &order, vector<ll> &rho, vector<bool> &inBase, vector<Field> &p)
{
	K.resize(K.size() + 1);
	P.resize(P.size() + 1);
	rho.resize(rho.size() + 1, -1);
	order.resize(order.size() + 1, INF);
	inBase.resize(inBase.size() + 1);
	p.resize(p.size() + 1);

	Q.X.resize(Q.row + 1);
	Q.row++;
	for (ll i = 0; i < Q.row; i++)
	{
		Q.X[i].resize(Q.col + 1);
	}
	Q.col++;

	vertex_number++;
	return vertex_number - 1;
}

//remove vertex from V
void remove_vertex(ll u, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, Matrix &C)
{
	/*
	cout << "-----remove vertex start-----" << endl;
	/////
	cout << "remove " << u << endl;
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;
	cout << "Ba_sub:";
	for (ll i = 0; i < Ba_sub.size(); i++)
	{
		cout << Ba_sub[i] << " ";
	}
	cout << endl;
	cout << "NBa_sub:";
	for (ll i = 0; i < NBa_sub.size(); i++)
	{
		cout << NBa_sub[i] << " ";
	}
	cout << endl;
	*/
	/////
	ll u_idx;
	if (inBase[u])
	{
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == u)
			{
				u_idx = i;
				break;
			}
		}
		Ba.erase(Ba.begin() + u_idx);
		Ba.shrink_to_fit();
		//cout << u_idx << endl; /////
		ll u_sub_idx;
		if (Ba_sub == Ba)
		{
			for (ll i = 0; i < Ba_sub.size(); i++)
			{
				if (Ba_sub[i] == u)
				{
					u_sub_idx = i;
					break;
				}
			}
			Ba_sub.erase(Ba_sub.begin() + u_sub_idx);
			Ba_sub.shrink_to_fit();
		}

		C.X.erase(C.X.begin() + u_idx);
		C.X.shrink_to_fit();
		C.row--;
	}
	else
	{
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == u)
			{
				u_idx = i;
				break;
			}
		}
		//cout << u_idx << endl; /////

		NBa.erase(NBa.begin() + u_idx);
		NBa.shrink_to_fit();

		if (NBa_sub == NBa)
		{
			ll u_sub_idx;
			for (ll i = 0; i < NBa_sub.size(); i++)
			{
				if (NBa_sub[i] == u)
				{
					u_sub_idx = i;
					break;
				}
			}
			NBa_sub.erase(NBa_sub.begin() + u_sub_idx);
			NBa_sub.shrink_to_fit();
		}

		for (ll i = 0; i < C.row; i++)
		{
			C.X[i].erase(C.X[i].begin() + u_idx);
			C.X[i].shrink_to_fit();
		}
		C.col--;
	}

	/*
	cout << "remove " << u << endl;
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;
	cout << "Ba_sub:";
	for (ll i = 0; i < Ba_sub.size(); i++)
	{
		cout << Ba_sub[i] << " ";
	}
	cout << endl;
	cout << "NBa_sub:";
	for (ll i = 0; i < NBa_sub.size(); i++)
	{
		cout << NBa_sub[i] << " ";
	}
	cout << endl;

	cout << "-----remove vertex end------" << endl;
	*/
}

//expand a blossom H
void Expand(Matrix &C, node *H, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, Tree &Bl)
{
	//cout << "-----Expand start------" << endl;
	if (!H->normal)
	{
		tree_delete(Bl, H->key);
	}
	else
	{
		ll t = H->tip;
		ll b = H->bud;

		/////
		//cout<<t<<" "<<b<<endl;
		/////

		pr p;
		ll t_idx, b_idx;
		if (inBase[t])
		{
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == t)
				{
					t_idx = i;
					break;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == b)
				{
					b_idx = i;
					break;
				}
			}

			p = pr(t_idx, b_idx);
		}
		else
		{
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == b)
				{
					b_idx = i;
					break;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == t)
				{
					t_idx = i;
					break;
				}
			}
			p = pr(b_idx, t_idx);
			//p.first=b; p.second=t;
		}

		Pivoting_around_p(C, p, Ba, NBa, Ba_sub, NBa_sub, inBase);

		//↓Search in Blossomの時だけ関係

		vector<node *> H_ancestors = ancestors(Bl, H);
		for (ll i = 0; i < H_ancestors.size(); i++)
		{
			//cout << "sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
			node *Hi = H_ancestors[i];
			ll idx;
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				ll k = Hi->ordering[j];
				if (k == t)
				{
					idx = j;
					break;
				}
			}
			auto itr = Hi->ordering.erase(Hi->ordering.begin() + idx);

			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				ll k = Hi->ordering[j];
				if (k == b)
				{
					idx = j;
					break;
				}
			}
			itr = Hi->ordering.erase(Hi->ordering.begin() + idx);
		}

		tree_delete(Bl, H->key);
		remove_vertex(t, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
		remove_vertex(b, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
	}
	//cout << "-----Expand end-----" << endl;
}

//judge if there exists an edge between u and v
bool existsEdge(ll u, ll v, Matrix &C, Matrix &Q, vector<ll> &Ba, vector<ll> &NBa, vector<bool> &inBase, vector<Field> &p)
{
	ll u_idx;
	ll v_idx;
	if (inBase[u] == inBase[v])
	{
		return false;
	}
	if (inBase[u])
	{
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == u)
			{
				u_idx = i;
				break;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == v)
			{
				v_idx = i;
				break;
			}
		}
		if (!C.X[u_idx][v_idx].isZero())
		{
			if (p[v] - p[u] == Q.X[u][v])
			{
				return true;
			}
		}
	}
	else
	{
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == v)
			{
				v_idx = i;
				break;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == u)
			{
				u_idx = i;
				break;
			}
		}
		if (!C.X[v_idx][u_idx].isZero())
		{
			if (p[u] - p[v] == Q.X[v][u])
			{
				return true;
			}
		}
	}
	return false;
}

set<ll> GaussJordan2(Matrix &A)
{
	//cout << "GJ start" << endl;

	ll c = A.col;
	ll r = A.row;

	set<ll> ret;

	ll cnt = 0;
	for (ll j = 0; j < c; j++)
	{
		Field mx;
		mx.setZero();
		ll pivot = j;
		for (ll i = cnt; i < r; i++)
		{
			if (mx < A.X[i][j] || mx < -A.X[i][j])
			{
				mx = std::max(A.X[i][j], -A.X[i][j]);
				pivot = i;
			}
		}

		if (mx.isZero())
		{
			continue;
		}
		else
		{
			ret.insert(j);
		}
		A.row_swap(cnt, pivot);

		Field a = A.X[cnt][j].inverse();

		//multiply by a so that A[j][j]=1
		A.row_multiply(cnt, a);

		for (ll i = 0; i < r; i++)
		{
			if (i != cnt)
			{
				Field b = A.X[i][j];
				A.row_plus(i, cnt, -b);
			}
		}
		cnt++;
	}
	//cout << "GJ end" << endl;

	return ret;
}

vector<ll> Base_Greedy(Matrix &A, vector<Field> &p)
{
	Matrix B = A;
	vector<pair<Field, ll>> M;
	for (ll i = 0; i < A.col; i++)
	{
		M.push_back(make_pair(p[i], i));
	}
	sort(M.begin(), M.end());

	/* for (ll i = 0; i < A.col; i++)
	{
		M[i].first.output();
		cout << " " << M[i].second << endl;
	} */

	for (ll j = 0; j < A.col; j++)
	{
		for (ll i = 0; i < A.row; i++)
		{
			B.X[i][j] = A.X[i][M[j].second];
		}
	}

	set<ll> C = GaussJordan2(B);
	vector<ll> ret;
	for (auto itr = C.begin(); itr != C.end(); ++itr)
	{
		ret.push_back(M[*itr].second);
	}

	return ret;
}

/*
void Blossom(ll v, ll u, Matrix &C, Matrix &Q, vector<vector<ll>> &P, Tree &Bl, vector<ll> &K, vector<ll> &order, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<Field> &p, vector<ll> &rho, queue<ll> &que, ll N)
{
	//cout << "Blossom start" << endl;
	//-----step1-----
	vector<ll> path_u;
	vector<ll> path_v;

	path_u = Path(u, P);
	reverse(path_u.begin(), path_u.end());
	path_v = Path(v, P);
	reverse(path_v.begin(), path_v.end());

	set<ll> st_u;
	for (ll i = 0; i < path_u.size(); i++)
	{
		st_u.insert(path_u[i]);
	}
	set<ll> st_v;
	for (ll i = 0; i < path_v.size(); i++)
	{
		st_v.insert(path_v[i]);
	}


	ll c;
	ll d;
	ll c_idx;
	ll d_idx;

	for (ll i = path_v.size() - 1; i >= 0; i--)
	{
		node *H;
		ll cur = path_v[i];

		set<ll> K_c;
		if (K[cur] > 0)
		{
			H = tree_search(Bl, K[cur]).first;
			for (ll j = 0; j < H->ordering.size(); j++)
			{
				K_c.insert(H->ordering[j]);
			}
			if (H->normal)
			{
				K_c.insert(H->bud);
			}
		}
		else
		{
			K_c.insert(cur);
			K_c.insert(mate(cur));
		}

		int ope = 0;
		for (auto itr = K_c.begin(); itr != K_c.end(); itr++)
		{
			if (st_u.find(*itr) != st_u.end())
			{
				ope = 1;
				break;
			}
		}

		/////
		cout << "K(c):";
		for (auto itr = K_c.begin(); itr != K_c.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		/////

		if (ope == 1)
		{
			c = path_v[i];
			c_idx = i;
			for (ll j = path_u.size() - 1; j >= 0; j--)
			{
				if (K_c.find(path_u[j]) != K_c.end())
				{
					d = path_u[j];
					d_idx = j;
					break;
				}
			}
			break;
		}
	}

		cout << "c is " << c << endl;
	cout<<"d is "<<d<<endl;

	//ここまでOK

	ll r = -1;
	ll r_idx = -1;
	set<ll> children;
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		if (isSingle(path_v[i], K))
		{
			children.insert(path_v[i]); //自身だけで大丈夫かも
			children.insert(mate(path_v[i]));
		}
		else
		{
			children.insert(K[path_v[i]]);
		}
	}
	for (ll i = d_idx + 1; i < path_u.size(); i++)
	{
		if (isSingle(path_u[i], K))
		{
			children.insert(path_u[i]);
			children.insert(mate(path_u[i]));
		}
		else
		{
			children.insert(K[path_u[i]]);
		}
	}

	set<ll> H_elements;
	node *H;
	if (c == d)
	{
		r = c;
		r_idx = c_idx;
		Field q_H;
		q_H.setZero();
		blossom_insert(Bl, blossom_number, children, q_H);
		blossom_number++;

		H = tree_search(Bl, blossom_number - 1 + Bl.num).first;

		//rを求めたい
		node *x = H->fchild;
		while (x != NULL)
		{
			//葉じゃない(blossom)なら
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				H_elements.insert(x->ordering[i]);
			}
			if (!isLeaf(Bl, x))
			{
				if (x->normal)
				{
					H_elements.insert(x->bud);
				}
			}

			x = x->next;
		}
	}
	else
	{
		if (isSingle(c, K))
		{
			children.insert(c);
			children.insert(mate(c));
		}
		else{
			children.insert(K[c]);
		}

		/////OK
		cout << "child node:";
		for (auto itr = children.begin(); itr != children.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		/////
		Field q_H;
		q_H.setZero();
		blossom_insert(Bl, blossom_number, children, q_H);
		blossom_number++;

		H = tree_search(Bl, blossom_number - 1 + Bl.num).first;

		//rを求めたい
		node *x = H->fchild;
		while (x != NULL)
		{
			//葉じゃない(blossom)なら
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				H_elements.insert(x->ordering[i]);
			}
			if (!isLeaf(Bl, x))
			{
				if (x->normal)
				{
					H_elements.insert(x->bud);
				}
			}
			x = x->next;
		}
		for (ll i = c_idx - 1; i >= 0; i--)
		{
			if (H_elements.find(path_v[i]) == H_elements.end())
			{
				r = path_v[i];
				r_idx = i;
				break;
			}
		}

	}
	H->normal = false;

	//-----step2-----
	ll g = path_v[r_idx + 1];
	ll h;
	if (r_idx + 1 < path_u.size())
	{
		h = path_u[r_idx + 1];
	}
	else
	{
		h = g;
	}

	if (!SourceLine_set(N, H_elements, inBase))
	{
		//cout << "There is no source line. New blossom is a normal blossom." << endl;
		H->tip = new_vertex(K, P, C, Q, order, rho, inBase, p);
		H->bud = new_vertex(K, P, C, Q, order, rho, inBase, p);
		ll t = H->tip;
		ll b = H->bud;
		H->normal = true;

		if (inBase[r] == 1 && inBase[g] == 0)
		{
			Ba.push_back(H->bud);
			NBa.push_back(H->tip);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			Ba_sub.push_back(H->bud);
			NBa_sub.push_back(H->tip);

			inBase[H->bud] = true;
			inBase[H->tip] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				//q(H)=0とするからこれでよい
				Q.X[Ba[i]][H->tip] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->tip][Ba[i]] = Q.X[Ba[i]][Ba[i]];

				Q.X[Ba[i]][H->bud] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->bud][Ba[i]] = Q.X[Ba[i]][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][H->tip] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->tip][NBa[i]] = Q.X[NBa[i]][NBa[i]];

				Q.X[NBa[i]][H->bud] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->bud][NBa[i]] = Q.X[NBa[i]][NBa[i]];
			}

			ll b_idx = Ba.size() - 1;
			ll t_idx = NBa.size() - 1;
			ll r_idx;
			ll g_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == r)
				{
					r_idx = i;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == g)
				{
					g_idx = i;
				}
			}

			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(NBa[i]) != H_elements.end())
				{
					C.X[b_idx][i] = C.X[r_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (H_elements.find(Ba[i]) == H_elements.end())
				{
					C.X[i][t_idx] = C.X[i][g_idx];
				}
			}
			C.X[b_idx][t_idx] = C.X[r_idx][g_idx];

            //update Q
			Q.X[r][b] = Q.X[r][r];
			p[b] = p[r] + Q.X[r][b];
			p[t] = p[b];
		}
		else if (inBase[r] == 0 && inBase[g] == 1)
		{
			Ba.push_back(H->tip);
			NBa.push_back(H->bud);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			Ba_sub.push_back(H->tip);
			NBa_sub.push_back(H->bud);

			inBase[H->tip] = true;
			inBase[H->bud] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				//q(H)=0とするからこれでよい
				Q.X[Ba[i]][H->tip] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->tip][Ba[i]] = Q.X[Ba[i]][Ba[i]];

				Q.X[Ba[i]][H->bud] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->bud][Ba[i]] = Q.X[Ba[i]][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][H->tip] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->tip][NBa[i]] = Q.X[NBa[i]][NBa[i]];

				Q.X[NBa[i]][H->bud] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->bud][NBa[i]] = Q.X[NBa[i]][NBa[i]];
			}

			ll t_idx = Ba.size() - 1;
			ll b_idx = NBa.size() - 1;
			ll r_idx;
			ll g_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == g)
				{
					g_idx = i;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == r)
				{
					r_idx = i;
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(Ba[i]) != H_elements.end())
				{
					C.X[i][b_idx] = C.X[i][r_idx];
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				if (H_elements.find(NBa[i]) == H_elements.end())
				{
					C.X[t_idx][i] = C.X[g_idx][i];
				}
			}
			C.X[t_idx][b_idx] = C.X[g_idx][r_idx];

            //update Q
			Q.X[r][b] = Q.X[r][r];
			p[b] = p[r] - Q.X[r][b];
			p[t] = p[b];
		}

		Pivoting_around_p(C, pr(Ba.size() - 1, NBa.size() - 1), Ba, NBa, Ba_sub, NBa_sub, inBase);
		//-----step3-----

		//replacing path 
		if (!P[mate(g)].empty())
		{
			P[mate(g)][0] = t;
		}
		else
		{
			node *H_g = tree_search(Bl, K[g]).first;
			auto itr = P[H_g->bud].begin() + 1;
			itr = P[H_g->bud].insert(itr, H->bud);
			itr++;
			itr = P[H_g->bud].insert(itr, H->tip);
		}

		if (!P[mate(h)].empty())
		{
			P[mate(h)][0] = t;
		}
		else
		{
			node *H_h = tree_search(Bl, K[h]).first;
			auto itr = P[H_h->bud].begin() + 1;
			itr = P[H_h->bud].insert(itr, H->bud);
			itr++;
			itr = P[H_h->bud].insert(itr, H->tip);
		}

		cout << "replacing path end" << endl;

		//label t with P(t)=P(r)bt
		P[t].push_back(r);
		P[t].push_back(b);

		for (ll i = 0; i < order.size(); i++)
		{
			if (order[i] < INF && order[i] > order[r])
			{
				order[i]++;
			}
		}
		order[H->tip] = order[r] + 1;
		ordering_number++;

		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			if (rho[*itr] == r)
			{
				rho[*itr] = H->tip;
			}
		}
		rho[b] = r;
		H_elements.insert(H->tip);
		
	}

	//-----step4-----
	//cout << endl;
	//cout << "-----step4-----" << endl; 
	vector<pr> unlabeled;
	vector<pr> labeled;


	for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
	{
		if (order[*itr] < INF)
		{
			labeled.push_back(pr(order[*itr], *itr));
		}
	}

	//orderの値でソート(降順)
	sort(labeled.begin(), labeled.end(), comp);
	//labeledの頂点のみ先に＜Hを決める
	for (ll i = labeled.size() - 1; i >= 0; i--)
	{
		H->ordering.push_back(labeled[i].second);
	}

	node *x = H->fchild;
	queue<pr> que_g;
	queue<pr> que_h;

	//Step4の3，4番目の条件でlabelする頂点について
	//1,2番目の条件によりtiはすでにlabelされ，さらにRHi(x)もさだまっているから
	//P(x)=RHi(x)でよい
	//この後で必要なのはqueに入れる順番決め
	while (x != NULL)
	{
		if (x->label == -1 && st_u.find(x->tip) != st_u.end())
		{
			int ope_h = 0;
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				ll k = x->ordering[i];
				if (order[k] < INF || k == x->tip)
				{
					continue;
				}

				P[k] = x->routing[k];
				order[k] = ordering_number;
				//↑仮にorderをつけておく．あとで順番を考慮して付け直す

				if (r == c && r == d)
				{
					if (k == h)
					{
						que_h.push(pr(order[rho[k]], k));
						ope_h = 1;
					}
					else if (ope_h == 1)
					{
						que_h.push(pr(order[rho[k]], k));
					}
					else
					{
						unlabeled.push_back(pr(order[rho[k]], k));
					}
				}
				else
				{
					unlabeled.push_back(pr(order[rho[k]], k));
				}
			}
		}
		else if (x->label == -1 && st_v.find(x->tip) != st_v.end())
		{
			int ope_g = 0;
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				ll k = x->ordering[i];

				if (order[k] < INF || k == x->tip)
				{
					continue;
				}

				P[k] = x->routing[k];
				order[k] = ordering_number;

				if (r == c && r == d)
				{
					if (k == g)
					{
						que_g.push(pr(order[rho[k]], k));
						ope_g = 1;
					}
					else if (ope_g == 1)
					{
						que_g.push(pr(order[rho[k]], k));
					}
					else
					{
						unlabeled.push_back(pr(order[rho[k]], k));
					}
				}
				else
				{
					unlabeled.push_back(pr(order[rho[k]], k));
				}
			}
		}
		x = x->next;
	}

	//label unlabeled vertices with P(x)
	for (ll i = d_idx + 1; i < path_u.size(); i++)
	{
		ll k = path_u[i];
		if (order[k] == INF)
		{
			cout << "label " << k << endl;
			if (i + 2 >= path_u.size())
			{
				P[k].push_back(v);
			}
			else
			{
				P[k].push_back(path_u[i + 2]);
			}
			P[k].push_back(path_u[i + 1]);

			unlabeled.push_back(pr(order[rho[k]], k));
		}
	}
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		ll k = path_v[i];
		if (order[k] == INF)
		{
			if (i + 2 >= path_v.size())
			{
				if (u != r)
				{
					P[k].push_back(u);
				}
				else
				{
					P[k].push_back(H->tip);
				}
			}
			else
			{
				P[k].push_back(path_v[i + 2]);
			}
			P[k].push_back(path_v[i + 1]);

			unlabeled.push_back(pr(order[rho[k]], k));
		}
    }

	if (r == c && r == d)
	{
		if (!que_g.empty())
		{
			pr a = que_g.front();
			que_g.pop();
			unlabeled.push_back(a);
		}
		if (!que_h.empty())
		{
			pr b = que_h.front();
			que_h.pop();
			unlabeled.push_back(b);
		}

		while (!que_g.empty())
		{
			pr c = que_g.front();
			que_g.pop();
			unlabeled.push_back(c);
		}
		while (!que_h.empty())
		{
			pr c = que_h.front();
			que_h.pop();
			unlabeled.push_back(c);
		}
	}

	stable_sort(unlabeled.begin(), unlabeled.end(), comp);

	for (ll i = 0; i < unlabeled.size(); i++)
	{
		que.push(unlabeled[i].second);
		order[unlabeled[i].second] = ordering_number;
		ordering_number++;
		H->ordering.push_back(unlabeled[i].second);
	}

	//-----step5-----
	H->label = 1;

	H->q.setZero();

	H->routing.resize(vertex_number);
	if (!SourceLine_set(N, H_elements, inBase))
	{
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			if (*itr != H->tip)
			{
				H->routing[*itr] = P[*itr];
			}
			else
			{
				H->routing[H->tip] = {H->tip};
			}
			K[*itr] = H->key;
		}
	}
	else
	{
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			H->routing[*itr] = P[*itr];

			K[*itr] = H->key;
		}
	}
	K[H->bud] = H->key;

	cout << "Blossom end" << endl;

}
*/
void Blossom(ll v, ll u, Matrix &C, Matrix &Q, vector<vector<ll>> &P, Tree &Bl, Tree &Bl2, vector<ll> &K, vector<ll> &order, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<Field> &p, vector<ll> &rho, queue<ll> &que, ll N)
{
	if (!overflow)
	{
		Blossom_num++;
	}

	//cout<<"---------------------"<<endl;
	//cout<<"Blossom("<<u<<","<<v<<")"<<endl;
	/*
	cout << "Blossom start" << endl;
	cout << "ordering_number:" << ordering_number << endl;
	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;

	for (ll i = 0; i < P.size(); i++)
	{
		cout << "P[" << i << "]:";
		for (ll j = 0; j < P[i].size(); j++)
		{
			cout << P[i][j] << " ";
		}
		cout << endl;
	}
	*/

	/* for (ll i = 0; i < K.size(); i++)
		{
			cout << "i:" << i << " K[i]:" << K[i] << endl;
		} */

	/* cout << "C:" << endl;
	C.output_matrix(); */
	//-----step1-----
	vector<ll> path_u;
	vector<ll> path_v;

	path_u = Path(u, P);
	reverse(path_u.begin(), path_u.end());

	/////
	/*
	cout << "path_u:";
	for (ll i = 0; i < path_u.size(); i++)
	{
		cout << path_u[i] << " ";
	}
	*/
	/////
	path_v = Path(v, P);
	reverse(path_v.begin(), path_v.end());

	/////
	/*
	cout << endl;
	cout << "path_v:";
	for (ll i = 0; i < path_v.size(); i++)
	{
		cout << path_v[i] << " ";
	}
	cout << endl;
	*/
	/////

	set<ll> st_u;
	for (ll i = 0; i < path_u.size(); i++)
	{
		st_u.insert(path_u[i]);
	}
	set<ll> st_v;
	for (ll i = 0; i < path_v.size(); i++)
	{
		st_v.insert(path_v[i]);
	}

	//cout<<"path_u.size:"<<path_u.size()<<endl;

	ll c;
	ll d;
	ll c_idx;
	ll d_idx;

	for (ll i = path_v.size() - 1; i >= 0; i--)
	{
		node *H;
		ll cur = path_v[i];
		//cout << "cur is " << cur << endl;
		//cout << "K[cur] is " << K[cur] << endl;

		set<ll> K_c;

		if (K[cur] > 0)
		{
			H = tree_search(Bl, K[cur]).first;

			//cout << "H->key is " << H->key << endl;

			/*
			for (ll j = 0; j < H->ordering.size(); j++)
			{
				//P(u)の頂点に含まれてたら
				if (st_u.find(H->ordering[j]) != st_u.end())
				{
					ope = 1;
					break;
				}
			}
			if(H->normal){
				if(st_u.find(H->bud)!=st_u.end()){
					ope=1;break;
				}
			}
			*/
			for (ll j = 0; j < H->ordering.size(); j++)
			{
				K_c.insert(H->ordering[j]);
			}
			if (H->normal)
			{
				K_c.insert(H->bud);
			}
		}
		else
		{
			//このときcurはsingle vertex
			//H = tree_search(Bl, -K[cur]).first;

			//cout << "H->key is " << H->key << endl;
			/*
			if(st_u.find(H->ordering[cur])!=st_u.end()||st_u.find(H->ordering[mate(cur)])!=st_u.end()){
				ope=1; break;
			}
			*/
			K_c.insert(cur);
			K_c.insert(mate(cur));
		}

		int ope = 0;
		for (auto itr = K_c.begin(); itr != K_c.end(); itr++)
		{
			if (st_u.find(*itr) != st_u.end())
			{
				ope = 1;
				break;
			}
		}

		/////
		/*
		cout << "K(c):";
		for (auto itr = K_c.begin(); itr != K_c.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		*/
		/////

		if (ope == 1)
		{
			c = path_v[i];
			//cout << "c is " << c << endl;
			c_idx = i;
			/*
			set<ll> st2;
			for (ll j = 0; j < H->ordering.size(); j++)
			{
				st2.insert(H->ordering[j]);
			}
			if(H->normal){
				st2.insert(H->bud);
			}
			
			for (ll j = path_u.size(); j >= 0; j--)
			{
				if (st_2.find(path_u[j]) != st2.end())
				{
					d = path_u[j];
					d_idx = j;
					break;
				}
			}
			*/
			/////
			/*
			cout<<"path_u:";
			for(ll j=path_u.size()-1;j>=0;j--)
			{
				cout<<path_u[j]<<" ";
			}
			cout<<endl;
			*/
			for (ll j = path_u.size() - 1; j >= 0; j--)
			{
				//cout<<"path_u[j]:"<<path_u[j]<<endl;
				if (K_c.find(path_u[j]) != K_c.end())
				{
					d = path_u[j];
					d_idx = j;
					break;
				}
			}
			//cout << "d is " << d << endl;
			break;
		}
	}
	//cout<<"c is "<<c<<endl;
	//cout<<"d is "<<d<<endl;

	//ここまでOK

	ll r = -1;
	ll r_idx = -1;
	set<ll> children;
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		if (isSingle(path_v[i], K))
		{
			children.insert(path_v[i]); //自身だけで大丈夫かも
			children.insert(mate(path_v[i]));
		}
		else
		{
			children.insert(K[path_v[i]]);
		}
	}
	for (ll i = d_idx + 1; i < path_u.size(); i++)
	{
		if (isSingle(path_u[i], K))
		{
			children.insert(path_u[i]);
			children.insert(mate(path_u[i]));
		}
		else
		{
			children.insert(K[path_u[i]]);
		}
	}

	set<ll> H_elements;
	node *H;
	if (c == d)
	{
		/////OK
		/*
		cout << "child node:";
		for (auto itr = children.begin(); itr != children.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		*/
		/////

		r = c;
		r_idx = c_idx;
		Field q_H;
		q_H.setZero();
		blossom_insert(Bl, blossom_number, children, q_H);
		blossom_number++;

		H = tree_search(Bl, blossom_number - 1 + Bl.num).first;

		//rを求めたい
		node *x = H->fchild;
		//set<ll> H_elements;
		while (x != NULL)
		{
			//葉じゃない(blossom)なら
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				H_elements.insert(x->ordering[i]);
			}
			if (!isLeaf(Bl, x))
			{
				if (x->normal)
				{
					H_elements.insert(x->bud);
				}
			}

			x = x->next;
		}

		/////
		/*
		cout << "key of new blossom is ";
		cout << H->key << endl;

		cout << "elements of new blossom: ";
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		*/
		/////
	}
	else
	{
		if (isSingle(c, K))
		{
			children.insert(c);
			children.insert(mate(c));
		}
		else
		{
			children.insert(K[c]);
		}

		/////OK
		/*
		cout << "child node:";
		for (auto itr = children.begin(); itr != children.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		*/
		/////

		Field q_H;
		q_H.setZero();
		blossom_insert(Bl, blossom_number, children, q_H);
		blossom_number++;

		H = tree_search(Bl, blossom_number - 1 + Bl.num).first;

		//rを求めたい
		node *x = H->fchild;
		//set<ll> H_elements;
		while (x != NULL)
		{
			//cout << endl;
			//cout << "x->key is " << x->key << endl;
			//葉じゃない(blossom)なら
			for (ll i = 0; i < x->ordering.size(); i++)
			{
				//cout << "x->ordering" << x->ordering[i] << " ";
				H_elements.insert(x->ordering[i]);
			}
			if (!isLeaf(Bl, x))
			{
				if (x->normal)
				{
					H_elements.insert(x->bud);
				}
			}
			x = x->next;
		}
		for (ll i = c_idx - 1; i >= 0; i--)
		{
			if (H_elements.find(path_v[i]) == H_elements.end())
			{
				r = path_v[i];
				r_idx = i;
				break;
			}
		}

		/////
		/*
		cout << "r is " << r << endl;
		cout << "key of new blossom is ";
		cout << H->key << endl;

		cout << "elements of new blossom: ";
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			cout << *itr << " ";
		}
		cout << endl;
		*/
		/////
	}
	H->normal = false;
	H->routing.resize(vertex_number + 2);

	//-----step2-----
	ll g = path_v[r_idx + 1];

	/////!!!!!!!!!hが存在しない可能性
	ll h;
	if (r_idx + 1 < path_u.size())
	{
		h = path_u[r_idx + 1];
	}
	else
	{
		h = g;
	}

	/////
	/*
	cout << "r is " << r << endl;
	cout << "g is " << g << endl;
	cout << "h is " << h << endl;
	*/
	/////

	if (!SourceLine_set(N, H_elements, inBase))
	{
		//cout << "There is no source line. New blossom is a normal blossom." << endl;
		//ll g = path_v[r_idx + 1];
		H->tip = new_vertex(K, P, C, Q, order, rho, inBase, p);
		H->bud = new_vertex(K, P, C, Q, order, rho, inBase, p);

		/////
		//↓Search _in_blossomの時だけ関係
		vector<node *> H_ancestors = ancestors(Bl2, H);
		for (ll i = 0; i < H_ancestors.size(); i++)
		{
			node *Hi = H_ancestors[i];
			Hi->ordering.push_back(H->tip);
			Hi->ordering.push_back(H->bud);
		}
		/////
		/*
		cout << "the tip of H is " << H->tip << " and the bud is " << H->bud << endl;
		*/
		/////

		ll t = H->tip;
		ll b = H->bud;
		H->normal = true;

		for (ll i = 0; i < Ba.size(); i++)
		{
			//q(H)=0とするからこれでよい
			Q.X[Ba[i]][H->tip] = Q.X[Ba[i]][Ba[i]];
			Q.X[H->tip][Ba[i]] = Q.X[Ba[i]][Ba[i]];

			Q.X[Ba[i]][H->bud] = Q.X[Ba[i]][Ba[i]];
			Q.X[H->bud][Ba[i]] = Q.X[Ba[i]][Ba[i]];
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			Q.X[NBa[i]][H->tip] = Q.X[NBa[i]][NBa[i]];
			Q.X[H->tip][NBa[i]] = Q.X[NBa[i]][NBa[i]];

			Q.X[NBa[i]][H->bud] = Q.X[NBa[i]][NBa[i]];
			Q.X[H->bud][NBa[i]] = Q.X[NBa[i]][NBa[i]];
		}

		//↓Search_in_Blossomの時だけ関係
		if (Bl.root->key != -1)
		{
			node *Hj = Bl.root;
			//cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
			ll tj = Hj->tip;
			//cout << Hj->key << endl;
			//Hj->q.output();
			//cout << endl;

			for (ll i = 0; i < Ba.size(); i++)
			{
				ll k = Ba[i];
				Q.X[H->tip][k] = Q.X[tj][k];
				Q.X[k][H->tip] = Q.X[k][tj];
				Q.X[H->bud][k] = Q.X[tj][k];
				Q.X[k][H->bud] = Q.X[k][tj];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				ll k = NBa[i];
				Q.X[H->tip][k] = Q.X[tj][k];
				Q.X[k][H->tip] = Q.X[k][tj];
				Q.X[H->bud][k] = Q.X[tj][k];
				Q.X[k][H->bud] = Q.X[k][tj];
			}
			Q.X[H->tip][tj].setZero();
			Q.X[tj][H->tip].setZero();
			Q.X[H->bud][tj].setZero();
			Q.X[tj][H->bud].setZero();
			Q.X[H->tip][H->tip] = Q.X[tj][tj];
			Q.X[H->bud][H->bud] = Q.X[tj][tj];
			Q.X[H->tip][H->bud].setZero();
			Q.X[H->bud][H->tip].setZero();
			/*	
				for(ll i=0;i<Ba.size();i++){
					ll k=Ba[i];
					Q.X[k][H->tip]+=Hj->q;
					Q.X[H->tip][k]+=Hj->q;
					Q.X[k][H->bud]+=Hj->q;
					Q.X[H->bud][k]+=Hj->q;
				}
				for(ll i=0;i<NBa.size();i++){
					ll k=NBa[i];
					Q.X[k][H->tip]+=Hj->q;
					Q.X[H->tip][k]+=Hj->q;
					Q.X[k][H->bud]+=Hj->q;
					Q.X[H->bud][k]+=Hj->q;
				}
				for(ll i=0;i<Hj->ordering.size();i++){
					ll k=Hj->ordering[i];
					Q.X[k][H->tip]-=Hj->q.Double();
					Q.X[H->tip][k]-=Hj->q.Double();
					Q.X[k][H->bud]-=Hj->q.Double();
					Q.X[H->bud][k]-=Hj->q.Double();
				}
				Q.X[H->bud][H->bud]=Hj->q;
				Q.X[H->tip][H->tip]=Hj->q;
				*/
			//cout << "^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		}

		if (inBase[r] == 1 && inBase[g] == 0)
		{
			Ba.push_back(H->bud);
			NBa.push_back(H->tip);

			//cout << "Ba:";
			/*
			for (ll i = 0; i < Ba.size(); i++)
			{
				cout << Ba[i] << " ";
			}
			cout << endl;
			cout << "Ba size is ";
			cout << Ba.size() << endl;
			*/

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;
			/////
			/*
			cout << "size of C update:" << endl;
			C.output_matrix();
			cout << endl;
			*/
			/////

			Ba_sub.push_back(H->bud);
			NBa_sub.push_back(H->tip);

			inBase[H->bud] = true;
			inBase[H->tip] = false;
			/*
			
			for (ll i = 0; i < Ba.size(); i++)
			{
				//q(H)=0とするからこれでよい
				Q.X[Ba[i]][H->tip] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->tip][Ba[i]] = Q.X[Ba[i]][Ba[i]];

				Q.X[Ba[i]][H->bud] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->bud][Ba[i]] = Q.X[Ba[i]][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][H->tip] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->tip][NBa[i]] = Q.X[NBa[i]][NBa[i]];

				Q.X[NBa[i]][H->bud] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->bud][NBa[i]] = Q.X[NBa[i]][NBa[i]];
			}
			*/

			/////
			/* cout << "Q was updated: " << endl;
			Q.output_matrix();
			cout << endl; */
			/////

			//Q.X[H->tip][H->tip].setZero();

			//update C
			/*
			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;
			*/

			ll b_idx = Ba.size() - 1;
			ll t_idx = NBa.size() - 1;
			ll r_idx;
			ll g_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == r)
				{
					r_idx = i;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == g)
				{
					g_idx = i;
				}
			}

			//多分このままでSearch_in_blossomでもいける
			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(NBa[i]) != H_elements.end())
				{
					C.X[b_idx][i] = C.X[r_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (H_elements.find(Ba[i]) == H_elements.end())
				{
					C.X[i][t_idx] = C.X[i][g_idx];
				}
			}
			C.X[b_idx][t_idx] = C.X[r_idx][g_idx];

			/* cout << "C was update: " << endl;
			C.output_matrix();
			cout << endl; */

			/////ここはOK

			/*
			for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
			{
				if (inBase[*itr] == 0)
				{
					ll idx = NBa_inverse[*itr];
					C.X[Ba.size() - 1][idx] = C.X[Ba_inverse[r]][idx];
				}
			}
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (H_elements.find(Ba[i]) != H_elements.end())
				{
					C.X[i][NBa.size() - 1] = C.X[i][NBa_inverse[g]];
				}
			}
			C.X[Ba.size() - 1][NBa.size() - 1] = C.X[Ba_inverse[r]][NBa_inverse[g]];
			*/

			//p.resize(p.size() + 2);

			//!!!!!Qの更新
			p[b] = p[r] + Q.X[r][b];
			p[t] = p[b];
		}
		else if (inBase[r] == 0 && inBase[g] == 1)
		{
			Ba.push_back(H->tip);
			NBa.push_back(H->bud);

			/*
			cout << "Ba:";
			for (ll i = 0; i < Ba.size(); i++)
			{
				cout << Ba[i] << " ";
			}
			cout << endl;
			/* cout << "Ba size is ";
			cout << Ba.size() << endl; */

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;
			/*
			cout << "size of C update:" << endl;
			C.output_matrix();
			cout << endl;
			*/

			Ba_sub.push_back(H->tip);
			NBa_sub.push_back(H->bud);

			inBase[H->tip] = true;
			inBase[H->bud] = false;

			/*
			for (ll i = 0; i < Ba.size(); i++)
			{
				//q(H)=0とするからこれでよい
				Q.X[Ba[i]][H->tip] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->tip][Ba[i]] = Q.X[Ba[i]][Ba[i]];

				Q.X[Ba[i]][H->bud] = Q.X[Ba[i]][Ba[i]];
				Q.X[H->bud][Ba[i]] = Q.X[Ba[i]][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][H->tip] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->tip][NBa[i]] = Q.X[NBa[i]][NBa[i]];

				Q.X[NBa[i]][H->bud] = Q.X[NBa[i]][NBa[i]];
				Q.X[H->bud][NBa[i]] = Q.X[NBa[i]][NBa[i]];
			}
			*/

			/////
			/* cout << "Q was updated: " << endl;
			Q.output_matrix();
			cout << endl; */
			/////

			//update C
			/*
			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;
			*/

			/*
			for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
			{
				if (inBase[*itr] == 1)
				{
					ll idx = Ba_inverse[*itr];
					C.X[idx][NBa.size() - 1] = C.X[idx][NBa_inverse[g]];
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (H_elements.find(NBa[i]) != H_elements.end())
				{
					C.X[Ba.size() - 1][i] = C.X[Ba_inverse[g]][i];
				}
			}
			C.X[Ba.size() - 1][NBa.size() - 1] = C.X[Ba_inverse[g]][NBa_inverse[r]];
			*/
			ll t_idx = Ba.size() - 1;
			ll b_idx = NBa.size() - 1;
			ll r_idx;
			ll g_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == g)
				{
					g_idx = i;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == r)
				{
					r_idx = i;
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(Ba[i]) != H_elements.end())
				{
					C.X[i][b_idx] = C.X[i][r_idx];
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				if (H_elements.find(NBa[i]) == H_elements.end())
				{
					C.X[t_idx][i] = C.X[g_idx][i];
				}
			}
			C.X[t_idx][b_idx] = C.X[g_idx][r_idx];

			/////
			/* cout << "C was update: " << endl;
			C.output_matrix();
			cout << endl; */
			/////

			//!!!!!Qの更新
			p[b] = p[r] - Q.X[r][b];
			p[t] = p[b];
		}

		/////
		/////

		Pivoting_around_p(C, pr(Ba.size() - 1, NBa.size() - 1), Ba, NBa, Ba_sub, NBa_sub, inBase);

		/*
		cout << "Ba:";
		for (ll i = 0; i < Ba.size(); i++)
		{
			cout << Ba[i] << " ";
		}
		cout << endl;
		cout << "NBa:";
		for (ll i = 0; i < NBa.size(); i++)
		{
			cout << NBa[i] << " ";
		}
		cout << endl;
		cout << "Ba_sub:";
		for (ll i = 0; i < Ba_sub.size(); i++)
		{
			cout << Ba_sub[i] << " ";
		}
		cout << endl;
		cout << "NBa_sub:";
		for (ll i = 0; i < NBa_sub.size(); i++)
		{
			cout << NBa_sub[i] << " ";
		}
		cout << endl;
		*/

		//-----step3-----
		//ll h = path_u[r_idx + 1];

		//gとｈはunlabeled
		/*
		auto itr = P[mate(g)].begin() + 1;
		P[mate(g)].insert(itr, H->bud);
		itr++;
		P[mate(g)].insert(itr, H->tip);

		if (h != g)
		{
			auto itr = P[mate(h)].begin() + 1;
			P[mate(h)].insert(itr, H->bud);
			itr++;
			P[mate(h)].insert(itr, H->tip);
		}
		*/

		//cout << "replacing path start";
		/* for (ll i = 0; i < K.size(); i++)
		{
			cout << "i:" << i << " K[i]:" << K[i] << endl;
		} */
		if (!P[mate(g)].empty())
		{
			P[mate(g)][0] = t;
		}
		else
		{
			vector<node *> des = descendant(Bl, H);

			/*
			node *H_g = tree_search(Bl, K[g]).first;
			auto itr = P[H_g->bud].begin() + 1;
			itr = P[H_g->bud].insert(itr, H->bud);
			itr++;
			itr = P[H_g->bud].insert(itr, H->tip);

			H->routing[H_g->bud].reserve(P[H_g->bud].size() - 2);
			itr = H->routing[H_g->bud].insert(H->routing[H_g->bud].begin(), P[H_g->bud].begin() + 2, P[H_g->bud].end());
			*/
		}

		/* for (ll i = 0; i < K.size(); i++)
		{
			cout << "i:" << i << " K[i]:" << K[i] << endl;
		} */
		if (!P[mate(h)].empty())
		{
			P[mate(h)][0] = t;
		}
		else
		{
			/*
			node *H_h = tree_search(Bl, K[h]).first;
			auto itr = P[H_h->bud].begin() + 1;
			itr = P[H_h->bud].insert(itr, H->bud);
			itr++;
			itr = P[H_h->bud].insert(itr, H->tip);

			H->routing[H_h->bud].reserve(P[H_h->bud].size() - 2);
			itr = H->routing[H_h->bud].insert(H->routing[H_h->bud].begin(), P[H_h->bud].begin() + 2, P[H_h->bud].end());
			*/
		}
		vector<node *> des = descendant(Bl, H);

		for (ll i = 0; i < des.size(); i++)
		{
			node *Hi = des[i];
			if (Hi->normal && Hi->label == -1)
			{
				if (P[Hi->bud][0] == r)
				{
					auto itr = P[Hi->bud].begin() + 1;
					itr = P[Hi->bud].insert(itr, H->bud);
					itr++;
					itr = P[Hi->bud].insert(itr, H->tip);

					H->routing[Hi->bud].reserve(P[Hi->bud].size() - 2);
					itr = H->routing[Hi->bud].insert(H->routing[Hi->bud].begin(), P[Hi->bud].begin() + 2, P[Hi->bud].end());
				}
			}
		}

		//cout << "replacing path end" << endl;

		//label t with P(t)=P(r)bt
		//P.resize(P.size() + 2);
		P[t].push_back(r);
		P[t].push_back(b);

		/////
		/* cout << endl;
		cout << "before step4" << endl;
		for (ll i = 0; i < P.size(); i++)
		{
			cout << "path i:";
			for (ll j = 0; j < P[i].size(); j++)
			{
				cout << P[i][j] << " ";
			}
			cout << endl;
		} */
		/////

		//extend the ordering

		//order.resize(order.size() + 2, INF);
		//order[H->bud] = INF;
		//cout << "before extending" << endl;
		/*
		cout << "ordeer:";
		for (ll i = 0; i < order.size(); i++)
		{
			cout << order[i] << " ";
		}
		cout << endl;
		*/
		/////
		for (ll i = 0; i < order.size(); i++)
		{
			if (order[i] < INF && order[i] > order[r])
			{
				order[i]++;
			}
		}
		order[H->tip] = order[r] + 1;
		ordering_number++;

		/////
		/*
		cout << endl;
		cout << "after extend the ordering <" << endl;
		cout << "order:";
		for (ll i = 0; i < order.size(); i++)
		{
			cout << order[i] << " ";
		}
		cout << endl;
		*/
		/////

		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			if (rho[*itr] == r)
			{
				rho[*itr] = H->tip;
			}
		}
		//rho[t]=-1;
		rho[b] = r;
		//rho.push_back(-1);
		//rho.push_back(r);

		/////
		/* cout << "rho:";
		for (ll i = 0; i < rho.size(); i++)
		{
			cout << rho[i] << " ";
		}
		cout << endl; */

		/* cout << "inBase:";
		for (ll i = 0; i < inBase.size(); i++)
		{
			if (inBase[i])
			{
				cout << "true"
					 << " ";
			}
			else
			{
				cout << "false"
					 << " ";
			}
		}
		cout << endl; */
		/////
		H_elements.insert(H->tip);
	}

	//-----step4-----
	//cout << endl;
	//cout << "-----step4-----" << endl;
	vector<pr> unlabeled;
	vector<pr> labeled;

	/////
	/*
	cout << "elements of H are:";
	for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
	{
		cout << *itr << " ";
	}
	cout << endl;
	*/
	/////

	for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
	{
		if (order[*itr] < INF)
		{
			labeled.push_back(pr(order[*itr], *itr));
		}
	}

	//orderの値でソート(降順)
	sort(labeled.begin(), labeled.end(), comp);
	//labeledの頂点のみ先に＜Hを決める
	for (ll i = labeled.size() - 1; i >= 0; i--)
	{
		H->ordering.push_back(labeled[i].second);
	}
	//cout << "labeled vertices" << endl;

	node *x = H->fchild;
	queue<pr> que_g;
	queue<pr> que_h;

	//Step4の3，4番目の条件でlabelする頂点について
	//1,2番目の条件によりtiはすでにlabelされ，さらにRHi(x)もさだまっているから
	//P(x)=RHi(x)でよい
	//この後で必要なのはqueに入れる順番決め
	/*
	vector<node*> blo=all_blossoms(Bl);
	for(ll i=0;i<blo.size();i++){
		if(!isLeaf(Bl,blo[i])){
			cout<<"key"<<blo[i]->key<<endl;
			cout<<"next key";
			if(blo[i]->next!=NULL){
				cout<<blo[i]->next->key<<endl;
			}
			else{
				cout<<"NULL"<<endl;
			}
		}
	}
	*/
	/////
	while (x != NULL)
	{
		//cout<<"a"<<endl;
		//cout<<x->key<<endl;
		if (x->label == 1)
		{
			if (x->normal && st_u.find(x->bud) != st_u.end())
			{

				ll b_idx;
				for (ll i = 0; i < path_u.size(); i++)
				{
					if (path_u[i] == x->bud)
					{
						b_idx = i;
						break;
					}
				}
				ll i = b_idx + 1;
				int ope = 0;
				if (order[x->bud] < INF)
				{
					continue;
				}
				while (i < path_u.size())
				{
					P[x->bud].push_back(path_u[i]);
					if (K[x->bud] != K[path_u[i]])
					{
						ope = 1;
						break;
					}
					i++;
				}
				if (ope == 0)
				{
					P[x->bud].push_back(v);
				}
				reverse(P[x->bud].begin(), P[x->bud].end());
				unlabeled.push_back(pr(order[rho[x->bud]], x->bud));
			}
			else if (x->normal && st_v.find(x->bud) != st_v.end())
			{
				ll b_idx;
				for (ll i = 0; i < path_v.size(); i++)
				{
					if (path_v[i] == x->bud)
					{
						b_idx = i;
						break;
					}
				}
				if (order[x->bud] < INF)
				{
					continue;
				}
				ll i = b_idx + 1;
				int ope = 0;
				while (i < path_v.size())
				{
					P[x->bud].push_back(path_v[i]);
					if (K[x->bud] != K[path_v[i]])
					{
						ope = 1;
						break;
					}
					i++;
				}
				if (ope == 0)
				{
					P[x->bud].push_back(u);
				}
				reverse(P[x->bud].begin(), P[x->bud].end());
				unlabeled.push_back(pr(order[rho[x->bud]], x->bud));
			}
		}
		if (x->label == -1 && st_u.find(x->tip) != st_u.end())
		{
			//vector<ll> vec;
			int ope_h = 0;
			for (ll I = 0; I < x->ordering.size(); I++)
			{
				ll k = x->ordering[I];
				if (order[k] < INF)
				{
					continue;
				}

				/*
				P[k].push_back(v);

				//+reverse(P(u|ti))
				for (ll j = path_u.size() - 1; path_u[j] != x->tip; j--)
				{
					P[k].push_back(path_u[j]);
				}
				//+RHi(x)
				vector<ll> R_k = Path(k, x->routing);
				reverse(R_k.begin(), R_k.end());
				P[k].reserve(P[k].size() + R_k.size());
				P[k].insert(P[k].end(), R_k.begin(), R_k.end());
				*/
				if (k == x->tip)
				{
					ll tx = x->tip;
					for (ll i = d_idx + 1; i < path_u.size(); i++)
					{
						if (path_u[i] == tx)
						{
							if (i + 2 >= path_u.size())
							{
								P[k].push_back(v);
							}
							else
							{
								P[k].push_back(path_u[i + 2]);
							}
							P[k].push_back(path_u[i + 1]);
							break;
						}
					}
				}
				else
				{
					P[k] = x->routing[k];
				}

				//order[k] = ordering_number;

				if (r == c && r == d)
				{
					/* 					
					if (k == g)
					{
						que_g.push(pr(order[rho[k]],k);
						ope_g=1;
					}
					else if (ope_g == 1)
					{
						que_g.push(pr(order[rho[k],k);
					}
					else
					{
						unlabeled.push_back(order[rho[k]], k);
					} */
					if (k == h)
					{
						que_h.push(pr(order[rho[k]], k));

						//cout << "que_h push:" << k << endl;
						ope_h = 1;
					}
					else if (ope_h == 1)
					{
						que_h.push(pr(order[rho[k]], k));
						//cout << "que_h push:" << k << endl;
					}
					else
					{
						unlabeled.push_back(pr(order[rho[k]], k));
					}
				}
				else
				{
					unlabeled.push_back(pr(order[rho[k]], k));
				}
			}
		}
		else if (x->label == -1 && st_v.find(x->tip) != st_v.end())
		{
			//vector<ll> vec;
			int ope_g = 0;
			for (ll I = 0; I < x->ordering.size(); I++)
			{
				ll k = x->ordering[I];

				if (order[k] < INF)
				{
					continue;
				}

				if (k == x->tip)
				{
					ll tx = x->tip;
					for (ll i = c_idx + 1; i < path_v.size(); i++)
					{
						if (path_v[i] == tx)
						{
							if (i + 2 >= path_v.size())
							{
								P[k].push_back(u);
							}
							else
							{
								P[k].push_back(path_v[i + 2]);
							}
							P[k].push_back(path_v[i + 1]);
							break;
						}
					}
				}
				else
				{
					P[k] = x->routing[k];
				}

				/*
					P[k].push_back(u);

					//+reverse(P(u|ti))
					for (ll j = path_v.size() - 1; path_v[j] != x->tip; j--)
					{
						P[k].push_back(path_u[j]);
					}
					//+RHi(x)
					vector<ll> R_k = Path(k, x->routing);
					reverse(R_k.begin(), R_k.end());
					P[k].insert(P[k].end(), R_k.begin(), R_k.end());
					*/

				//order[k] = ordering_number;

				if (r == c && r == d)
				{
					if (k == g)
					{
						que_g.push(pr(order[rho[k]], k));

						//cout << "que_g push:" << k << endl;
						ope_g = 1;
					}
					else if (ope_g == 1)
					{
						que_g.push(pr(order[rho[k]], k));
						//cout << "que_g push:" << k << endl;
					}
					else
					{
						unlabeled.push_back(pr(order[rho[k]], k));
					}
				}
				else
				{
					unlabeled.push_back(pr(order[rho[k]], k));
				}
			}
		}
		/*
		cout << "x key:" << x->key << endl;
		cout << "x next key:";
		if (x->next == NULL)
		{
			cout << "NULL" << endl;
		}
		else
		{
			cout << x->next->key << endl;
		}
		*/

		x = x->next;
	}

	//unlabeledの頂点をP(x)でラベル
	for (ll i = d_idx + 1; i < path_u.size(); i++)
	{
		ll k = path_u[i];
		if (order[k] == INF && isSingle(k, K))
		{
			/*
			P[k].push_back(v);
			for (ll j = path_u.size() - 1; j > i; j--)
			{
				P[k].push_back(path_u[j]);
			}
			*/
			//cout << "label " << k << endl;
			if (i + 2 >= path_u.size())
			{
				P[k].push_back(v);
			}
			else
			{
				P[k].push_back(path_u[i + 2]);
			}
			P[k].push_back(path_u[i + 1]);

			unlabeled.push_back(pr(order[rho[k]], k));
		}
	}
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		ll k = path_v[i];
		if (order[k] == INF && isSingle(k, K))
		{
			/*
			P[k].push_back(u);
			for (ll j = path_v.size() - 1; j > i; j--)
			{
				P[k].push_back(path_v[j]);
			}
			unlabeled.push_back(pr(order[rho[k]], k));
			*/
			//cout << "label " << k << endl;
			if (i + 2 >= path_v.size())
			{
				if (u != r)
				{
					P[k].push_back(u);
				}
				else
				{ //normal line 一本のみの花をつくるとき
					P[k].push_back(H->tip);
				}
			}
			else
			{
				P[k].push_back(path_v[i + 2]);
			}
			P[k].push_back(path_v[i + 1]);

			unlabeled.push_back(pr(order[rho[k]], k));
		}
	}

	/////OK
	/* cout << "after labeling vertices in new blossom H" << endl;
	/////
	for (ll i = 0; i < P.size(); i++)
	{
		cout << "path " << i << ":";
		for (ll j = 0; j < P[i].size(); j++)
		{
			cout << P[i][j] << " ";
		}
		cout << endl;
	} */
	/////
	/////

	if (r == c && r == d)
	{
		//cout << "except for que_g que_h:";

		/*
		for (ll i = 0; i < unlabeled.size(); i++)
		{
			cout << unlabeled[i].second << " ";
		}

		cout << endl;
		*/
		if (!que_g.empty())
		{
			pr a = que_g.front();
			que_g.pop();
			unlabeled.push_back(a);
		}
		if (!que_h.empty())
		{
			pr b = que_h.front();
			que_h.pop();
			unlabeled.push_back(b);
		}

		while (!que_g.empty())
		{
			pr c = que_g.front();
			que_g.pop();
			unlabeled.push_back(c);
		}
		while (!que_h.empty())
		{
			pr c = que_h.front();
			que_h.pop();
			unlabeled.push_back(c);
		}
	}

	/* cout << "vector unlabeled:" << endl;
	for (ll i = 0; i < unlabeled.size(); i++)
	{
		cout << "vertex:" << unlabeled[i].second << " order[rho[vertex]]:" << unlabeled[i].first << endl;
	} */

	stable_sort(unlabeled.begin(), unlabeled.end(), comp);

	/* cout << "vector unlabeled after sorting:" << endl;
	for (ll i = 0; i < unlabeled.size(); i++)
	{
		cout << "vertex:" << unlabeled[i].second << " order[rho[vertex]]:" << unlabeled[i].first << endl;
	} */

	for (ll i = 0; i < unlabeled.size(); i++)
	{
		que.push(unlabeled[i].second);
		order[unlabeled[i].second] = ordering_number;
		ordering_number++;
		H->ordering.push_back(unlabeled[i].second);
	}

	//-----step5-----
	H->label = 1;

	H->q.setZero();

	if (!SourceLine_set(N, H_elements, inBase))
	{
		/*
		H->routing.resize(vertex_number);
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			vector<ll> path_from_b = Path(*itr, P);
			for (ll i = 0; path_from_b[i] != H->bud; i++)
			{
				H->routing[*itr].push_back(path_from_b[i]);
			}
			reverse(H->routing[*itr].begin(), H->routing[*itr].end());

			//K(v)の更新
			K[*itr] = H->key;
		}
		*/
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			//if(P[*itr][0]==r){
			//	H->routing[*itr].resize(P[*itr].size()-2);
			//	H->routing[*itr].insert(H->routing[*itr].begin(),P[*itr].begin()+2,P[*itr].end());
			//}
			//else{
			if (H->routing[*itr].empty())
			{
				if (*itr != H->tip)
				{
					H->routing[*itr] = P[*itr];
				}
				else
				{
					H->routing[H->tip] = {H->tip};
				}

				//}
			}
			K[*itr] = H->key;
		}
		K[H->bud] = H->key;
	}
	else
	{
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
			H->routing[*itr] = P[*itr];

			K[*itr] = H->key;
		}
	}
	//K[H->bud] = H->key;

	//////
	/*
	cout<<"C[1697][1672]";
	if(vertex_number>=1697)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1697){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1672){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=1){C.X[idx_1][idx_2].output();}
	}
	cout<<endl;
	cout<<"C[1673][1255]";
	if(vertex_number>=1673)
	{
		ll idx_1=-1;ll idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1673){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){
			C.X[idx_1][idx_2].output();
		}
	}
	cout<<endl;
	cout<<"C[1697][1255]";
	if(vertex_number>=1697)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1697){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1)C.X[idx_1][idx_2].output();
	}
	cout<<endl;
	cout<<"C[1721][1255]";
	if(vertex_number>=1721)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1721){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1)C.X[idx_1][idx_2].output();
	}
	cout<<endl;
	cout<<"C[1723][1255]";
	if(vertex_number>=1723)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1723){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){
			C.X[idx_1][idx_2].output();
		}
	}
	cout<<endl;
	cout<<"---------------------"<<endl;
	*/
	//////

	/////
	/*
	cout << "K:";
	for (ll i = 0; i < vertex_number; i++)
	{
		cout << K[i] << " ";
	}
	cout << endl;
	cout << endl;
	cout << "order:";
	for (ll i = 0; i < vertex_number; i++)
	{
		cout << order[i] << " ";
	}
	cout << endl;
	cout << endl;

	cout << "H->key:" << H->key << endl;
	cout << "H->label:" << H->label << endl;
	cout << "H is ";
	if (H->normal)
	{
		cout << "normal" << endl;
	}
	else
	{
		cout << "source" << endl;
	}
	cout << "q(H) is ";
	H->q.output();
	cout << endl;
	cout << endl;

	cout << "ordering of H:";
	for (ll i = 0; i < H->ordering.size(); i++)
	{
		cout << H->ordering[i] << " ";
	}
	cout << endl;
	cout << endl;

	cout << "routing of H" << endl;
	for (ll i = 0; i < H->routing.size(); i++)
	{
		cout << "RH[" << i << "]:";
		for (ll j = 0; j < H->routing[i].size(); j++)
		{
			cout << H->routing[i][j] << " ";
		}
		cout << endl;
	}
	*/

	/////
	/*
	cout << "p was updated:" << endl;
	for (ll i = 0; i < p.size(); i++)
	{
		cout << i << " ";
		p[i].output();
		cout << endl;
	}

	cout << "Blossom end" << endl;
	*/

	//!!!!!update G
}

void Graft(ll v, node *Hi, Matrix &C, Matrix &Q, vector<vector<ll>> &P, Tree &Bl, Tree &Bl2, vector<ll> &K, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, vector<Field> &p, vector<ll> &rho, queue<ll> &que)
{
	if (!overflow)
	{
		Graft_num++;
	}
	//cout << "Graft start" << endl;
	//node *Hi = tree_search(Bl, k).first;
	ll bi = Hi->bud;
	ll ti = Hi->tip;

	ll r = v;
	ll g = bi;

	ll t = new_vertex(K, P, C, Q, order, rho, inBase, p);
	ll b = new_vertex(K, P, C, Q, order, rho, inBase, p);

	set<ll> H_elements;
	for (ll i = 0; i < Hi->ordering.size(); i++)
	{
		H_elements.insert(Hi->ordering[i]);
	}
	H_elements.insert(bi);
	H_elements.insert(t);

	for (ll i = 0; i < Ba.size(); i++)
	{
		//q(H)=epsとするからこれでよい
		Q.X[Ba[i]][t] = Q.X[Ba[i]][ti];
		Q.X[t][Ba[i]] = Q.X[ti][Ba[i]];
		//Hiはマキシマルブロッサムなのでｂはblossomに属さない
		Q.X[Ba[i]][b] = Q.X[Ba[i]][bi];
		Q.X[b][Ba[i]] = Q.X[bi][Ba[i]];
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		Q.X[NBa[i]][t] = Q.X[NBa[i]][ti];
		Q.X[t][NBa[i]] = Q.X[ti][NBa[i]];

		Q.X[NBa[i]][b] = Q.X[NBa[i]][bi];
		Q.X[b][NBa[i]] = Q.X[bi][NBa[i]];
	}
	Q.X[t][ti].setZero();
	Q.X[ti][t].setZero();
	Q.X[t][t] = Q.X[ti][ti];
	Q.X[t][b] = Q.X[ti][bi];
	Q.X[b][t] = Q.X[bi][ti];

	//Blossom step2と同様b,t導入
	if (inBase[r] == 1 && inBase[g] == 0)
	{
		Ba.push_back(b);
		NBa.push_back(t);

		C.X.resize(C.row + 1);
		C.row++;
		for (ll i = 0; i < C.row; i++)
		{
			C.X[i].resize(C.col + 1);
		}
		C.col++;

		Ba_sub.push_back(b);
		NBa_sub.push_back(t);

		inBase[b] = true;
		inBase[t] = false;

		/* 	for (ll i = 0; i < Ba.size(); i++)
		{
			//q(H)=epsとするからこれでよい
			Q.X[Ba[i]][t] = Q.X[Ba[i]][ti];
			Q.X[t][Ba[i]] = Q.X[ti][Ba[i]];
			//Hiはマキシマルブロッサムなのでｂはblossomに属さない
			Q.X[Ba[i]][b] = Q.X[Ba[i]][bi];
			Q.X[b][Ba[i]] = Q.X[bi][Ba[i]];
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			Q.X[NBa[i]][t] = Q.X[NBa[i]][ti];
			Q.X[t][NBa[i]] = Q.X[ti][NBa[i]];

			Q.X[NBa[i]][b] = Q.X[NBa[i]][bi];
			Q.X[b][NBa[i]] = Q.X[bi][NBa[i]];
		}
		Q.X[t][ti].setZero();
		Q.X[ti][t].setZero();
		Q.X[t][t] = Hi->q;
		Q.X[t][b] = Hi->q;
		Q.X[b][t] = Hi->q; */

		ll b_idx = Ba.size() - 1;
		ll t_idx = NBa.size() - 1;
		ll r_idx;
		ll g_idx;
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == r)
			{
				r_idx = i;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == g)
			{
				g_idx = i;
			}
		}
		for (ll i = 0; i < C.col; i++)
		{
			//if NBa[i] in H\B*
			if (H_elements.find(NBa[i]) != H_elements.end())
			{
				C.X[b_idx][i] = C.X[r_idx][i];
			}
		}
		for (ll i = 0; i < C.row; i++)
		{
			if (H_elements.find(Ba[i]) == H_elements.end())
			{
				C.X[i][t_idx] = C.X[i][g_idx];
			}
		}
		C.X[b_idx][t_idx] = C.X[r_idx][g_idx];

		//update p
		p[b] = p[r] + Q.X[r][b];
		p[t] = p[b];
	}
	else if (inBase[r] == 0 && inBase[g] == 1)
	{
		Ba.push_back(t);
		NBa.push_back(b);

		C.X.resize(C.row + 1);
		C.row++;
		for (ll i = 0; i < C.row; i++)
		{
			C.X[i].resize(C.col + 1);
		}
		C.col++;

		Ba_sub.push_back(t);
		NBa_sub.push_back(b);

		inBase[t] = true;
		inBase[b] = false;

		/* for (ll i = 0; i < Ba.size(); i++)
		{
			//q(H)=eps
			Q.X[Ba[i]][t] = Q.X[Ba[i]][ti];
			Q.X[t][Ba[i]] = Q.X[ti][Ba[i]];

			Q.X[Ba[i]][b] = Q.X[Ba[i]][Ba[i]];
			Q.X[b][Ba[i]] = Q.X[Ba[i]][Ba[i]];
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			Q.X[NBa[i]][t] = Q.X[NBa[i]][ti];
			Q.X[t][NBa[i]] = Q.X[ti][NBa[i]];

			Q.X[NBa[i]][b] = Q.X[NBa[i]][NBa[i]];
			Q.X[b][NBa[i]] = Q.X[NBa[i]][NBa[i]];
		}
		//Q[ti][ti]は本来０だがQ[i][i]を変えているため調整が必要
		//bに関してはHiがmaximalなので問題なし
		Q.X[t][ti].setZero();
		Q.X[ti][t].setZero();
		Q.X[t][t] = Hi->q;
		Q.X[b][t] = Hi->q;
		Q.X[ti][bi] = Hi->q; */

		ll t_idx = Ba.size() - 1;
		ll b_idx = NBa.size() - 1;
		ll r_idx;
		ll g_idx;
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == g)
			{
				g_idx = i;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == r)
			{
				r_idx = i;
			}
		}
		for (ll i = 0; i < C.row; i++)
		{
			//if NBa[i] in H\B*
			if (H_elements.find(Ba[i]) != H_elements.end())
			{
				C.X[i][b_idx] = C.X[i][r_idx];
			}
		}
		for (ll i = 0; i < C.col; i++)
		{
			if (H_elements.find(NBa[i]) == H_elements.end())
			{
				C.X[t_idx][i] = C.X[g_idx][i];
			}
		}
		C.X[t_idx][b_idx] = C.X[g_idx][r_idx];

		//update p
		p[b] = p[r] - Q.X[r][b];
		p[t] = p[b];
	}

	Pivoting_around_p(C, pr(Ba.size() - 1, NBa.size() - 1), Ba, NBa, Ba_sub, NBa_sub, inBase);

	//label t with P(t)=P(v)bt
	P[t].push_back(v);
	P[t].push_back(b);

	//extend the ordering of < so that t is just after v
	for (ll i = 0; i < order.size(); i++)
	{
		if (order[i] < INF && order[i] > order[v])
		{
			order[i]++;
		}
	}
	order[t] = order[v] + 1;
	ordering_number++;

	rho[b] = v;

	//-----step2-----
	//cout << "----------step2-----------" << endl;
	for (ll i = 0; i < Hi->ordering.size(); i++)
	{
		ll k = Hi->ordering[i];
		if (k == ti)
		{
			continue;
		}
		if (Hi->routing[k][0] == ti)
		{
			Hi->routing[k][0] = t;
		}
		P[k] = Hi->routing[k];
		order[k] = ordering_number;
		ordering_number++;

		que.push(k);
	}

	//-----step3-----
	//cout << "-----------step3---------" << endl;
	Hi->label = 1;
	Hi->ordering[0] = t;
	//!!!define RH(x) for each x in H\{bi}
	Hi->routing.resize(vertex_number);
	Hi->routing[t] = {t};

	Hi->bud = b;
	Hi->tip = t;
	K[b] = Hi->key;
	K[t] = Hi->key;

	//-----step4-----
	//cout << "-----------step4---------" << endl;
	Field eps = Hi->q;

	if (inBase[bi] == 1)
	{
		p[bi] += eps;
		p[t] += eps;
	}
	else
	{
		p[bi] -= eps;
		p[t] -= eps;
	}

	pr pivot;
	ll ti_idx, bi_idx;
	if (inBase[ti])
	{
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == ti)
			{
				ti_idx = i;
				break;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == bi)
			{
				bi_idx = i;
				break;
			}
		}
		pivot = pr(ti_idx, bi_idx);
	}
	else
	{
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == bi)
			{
				bi_idx = i;
				break;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == ti)
			{
				ti_idx = i;
				break;
			}
		}
		pivot = pr(bi_idx, ti_idx);
	}

	Pivoting_around_p(C, pivot, Ba, NBa, Ba_sub, NBa_sub, inBase);

	//↓Search_in_blossomの時だけ関係
	/////////
	vector<node *> Hi_ancestors = ancestors(Bl2, Hi);
	for (ll j = 0; j < Hi_ancestors.size(); j++)
	{
		node *Hj = Hi_ancestors[j];
		ll idx;
		for (ll k = 0; k < Hj->ordering.size(); k++)
		{
			if (Hj->ordering[k] == ti)
			{
				idx = k;
				break;
			}
		}
		auto itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
		for (ll k = 0; k < Hj->ordering.size(); k++)
		{
			if (Hj->ordering[k] == bi)
			{
				idx = k;
				break;
			}
		}
		itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
		Hj->ordering.push_back(b);
		Hj->ordering.push_back(t);
	}
	////////
	/*
	cout << "p was updated:" << endl;
	for (ll i = 0; i < p.size(); i++)
	{
		cout << i << " ";
		p[i].output();
		cout << endl;
	}
	*/

	remove_vertex(ti, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
	remove_vertex(bi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
}

vector<ll> Search(Matrix &C, vector<bool> &inBase, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, Matrix &Q, vector<ll> &K, Tree &Bl, Tree &Bl2, vector<Field> &p, vector<ll> &order, vector<ll> &rho, vector<vector<ll>> &P, ll N)
{
	if (!overflow)
		Search_num++;
	else
		return {};

	//cout<<"feasibility check"<<endl;
	/////
	for (ll i = 0; i < Ba.size(); i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			if (!C.X[i][j].isZero())
			{
				if (p[NBa[j]] - p[Ba[i]] < Q.X[Ba[i]][NBa[j]])
				{
					cout << "infeasible" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	vector<node *> L = all_blossoms(Bl);
	for (ll I = 0; I < L.size(); I++)
	{
		if (L[I]->normal)
		{
			ll t = L[I]->tip;
			ll b = L[I]->bud;

			ll t_idx, b_idx;
			if (inBase[t] && !inBase[b])
			{
				for (ll i = 0; i < Ba.size(); i++)
				{
					if (Ba[i] == t)
					{
						t_idx = i;
						break;
					}
				}
				for (ll i = 0; i < NBa.size(); i++)
				{
					if (NBa[i] == b)
					{
						b_idx = i;
						break;
					}
				}
				if (!C.X[t_idx][b_idx].isZero())
				{
					if (p[b] - p[t] != L[I]->q)
					{
						cout << "infeasible" << endl;
						exit(EXIT_FAILURE);
					}
				}
			}
			else if (inBase[b] && !inBase[t])
			{
				for (ll i = 0; i < Ba.size(); i++)
				{
					if (Ba[i] == b)
					{
						b_idx = i;
						break;
					}
				}
				for (ll i = 0; i < NBa.size(); i++)
				{
					if (NBa[i] == t)
					{
						t_idx = i;
						break;
					}
				}
				if (!C.X[b_idx][t_idx].isZero())
				{
					if (p[t] - p[b] != L[I]->q)
					{
						cout << "infeasible" << endl;
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}
	//cout<<"--------Search start-------------"<<endl;
	/////
	/*
	cout<<"C[1697][1672]";
	if(vertex_number>=1697)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1697){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1672){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){C.X[idx_1][idx_2].output();}
	}
	cout<<endl;
	cout<<"C[1673][1255]";
	if(vertex_number>=1673)
	{
		ll idx_1=-1;ll idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1673){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){
			C.X[idx_1][idx_2].output();
			cout<<endl;
		}
	}
	cout<<"C[1721][1255]";
	if(vertex_number>=1721)
	{
		ll idx_1=-1; ll idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1721){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){
			C.X[idx_1][idx_2].output();
			cout<<endl;
		}
	}
	cout<<"C[1723][1255]";
	if(vertex_number>=1723)
	{
		ll idx_1=-1,idx_2=-1;
		for(ll i=0;i<Ba.size();i++){
			if(Ba[i]==1723){
				idx_1=i;break;
			}
		}
		for(ll i=0;i<NBa.size();i++){
			if(NBa[i]==1255){
				idx_2=i;break;
			}
		}
		if(idx_1!=-1&&idx_2!=-1){
				C.X[idx_1][idx_2].output();
		cout<<endl;
		}
	
	}
	*/
	/////

	/////
	/*
	cout << endl;
	cout << "---------------Search start----------------------" << endl;
	C.output_matrix();

	cout << "Ba:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;
	cout << "Ba_sub:";
	for (ll i = 0; i < Ba_sub.size(); i++)
	{
		cout << Ba_sub[i] << " ";
	}
	cout << endl;
	cout << "NBa_sub:";
	for (ll i = 0; i < NBa_sub.size(); i++)
	{
		cout << NBa_sub[i] << " ";
	}
	cout << endl;
	*/

	//cout << "edges:" << endl;
	/*
	for (ll i = 0; i < Ba.size(); i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			ll v = Ba[i];
			ll u = NBa[j];

			ll v_idx, u_idx;
			for (ll k = 0; k < Ba.size(); k++)
			{
				if (Ba[k] == v)
				{
					v_idx = k;
					break;
				}
			}
			for (ll k = 0; k < NBa.size(); k++)
			{
				if (NBa[k] == u)
				{
					u_idx = k;
					break;
				}
			}
			if (!C.X[v_idx][u_idx].isZero())
			{
				if (p[u] - p[v] == Q.X[u][v])
				{
					cout << "tight:" << v << " " << u << endl;
				}
				else
				{
					cout << "no tight:" << v << " " << u << endl;
					p[v].output();
					cout << " ";
					p[u].output();
					cout << " ";
					Q.X[u][v].output();
					cout << endl;
				}
			}
		}
	}
	*/

	//cout << "Q:" << endl;
	//Q.output_matrix();
	//cout << endl;

	/* for (ll i = 0; i < K.size(); i++)
		{
			cout << "i:" << i << " K[i]:" << K[i] << endl;
		} */

	vector<ll> path;
	ll Num = inBase.size();

	queue<ll> que;

	//K[i] -i: single, number(0-|blossoms|) : maximal blossom
	P.clear();
	P.resize(Num);

	K.resize(Num);
	for (ll i = 0; i < K.size(); i++)
	{
		K[i] = -i;
	}

	rho.clear();
	rho.resize(Num, -1);

	order.clear();
	//initialize ordering with INF
	order.resize(Num, INF);
	ordering_number = 0;

	//initialization K
	blossom_initialize(Bl, K);

	/////
	/*
	if(Ba.size()!=Ba_sub.size())
	{
	for(ll i=0;i<Ba_sub.size();i++){
		cout<<Ba_sub[i]<<" "<<K[Ba_sub[i]]<<endl;
	}
	for(ll i=0;i<NBa_sub.size();i++){
		cout<<NBa_sub[i]<<" "<<K[NBa_sub[i]]<<endl;
	}
	}
	*/
	/////
	/*
	cout << "K:";
	for (ll i = 0; i < K.size(); i++)
	{
		cout << K[i] << " ";
	}
	cout << endl;
	*/
	/////

	/* for (ll i = 0; i < K.size(); i++)
		{
			cout << "i:" << i << " K[i]:" << K[i] << endl;
		} */

	set<pr> source_lines;

	for (ll i = 0; i < Ba_sub.size(); i++)
	{
		ll k = Ba_sub[i];
		if (!inBase[mate(k)] && isSingle(k, K))
		{
			P[k].push_back(k);
			order[k] = ordering_number;
			ordering_number++;
			que.push(k);
			pr sourceline;
			sourceline = pr(k, mate(k));
			source_lines.insert(sourceline);
		}
	}
	for (ll i = 0; i < NBa_sub.size(); i++)
	{
		ll k = NBa_sub[i];
		if (inBase[mate(k)] && isSingle(k, K))
		{
			P[k].push_back(k);
			order[k] = ordering_number;
			ordering_number++;
			que.push(k);
			pr sourceline;
			sourceline = pr(mate(k), k);
			source_lines.insert(sourceline);
		}
	}

	//unlabeled maximal source blossom
	vector<node *> MaximalBlossom = maximal_blossom(Bl);
	for (ll i = 0; i < MaximalBlossom.size(); i++)
	{

		node *x = MaximalBlossom[i];
		/////
		/*
		cout << "routing:" << endl;
		for (ll cnt = 0; cnt < x->ordering.size(); cnt++)
		{
			ll vertex = x->ordering[cnt];
			cout << "RH(" << vertex << "):";
			for (ll cnt2 = 0; cnt2 < x->routing[vertex].size(); cnt2++)
			{
				cout << x->routing[vertex][cnt2] << " ";
			}
			cout << endl;
		}
		*/
		/////
		if (!x->normal)
		{
			//labeling maximal blossom with plus
			x->label = 1;
			for (ll j = 0; j < x->ordering.size(); j++)
			{
				P[x->ordering[j]].reserve(x->routing[j].size());
				P[x->ordering[j]] = x->routing[x->ordering[j]];
				/* for (ll k = 0; k < x->routing[x->ordering[j]].size(); k++)
				{
					P[x->ordering[j]].push_back(x->routing[x->ordering[j]][k]);
				} */
				que.push(x->ordering[j]);
				order[x->ordering[j]] = ordering_number;
				ordering_number++;
			}
		}
	}

	//cout << "aa" << endl;
	for (auto itr = source_lines.begin(); itr != source_lines.end(); itr++)
	{
		pr sourceline = *itr;
		ll x = sourceline.first;
		ll y = sourceline.second;
		//cout << x << " " << y << endl;
		if (x >= N || y >= N)
		{
			continue;
		}
		ll x_idx;
		ll y_idx;
		for (ll i = 0; i < Ba.size(); i++)
		{
			if (Ba[i] == x)
			{
				x_idx = i;
				break;
			}
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			if (NBa[i] == y)
			{
				y_idx = i;
				break;
			}
		}
		if (!C.X[x_idx][y_idx].isZero())
		{
			//cout << "aaa" << endl;
			if (p[y] - p[x] == Q.X[x][y])
			{
				//cout << "aaaa" << endl;
				Blossom(x, y, C, Q, P, Bl, Bl2, K, order, Ba, NBa, Ba_sub, NBa_sub, inBase, p, rho, que, N);
			}
		}
	}
	/////
	//cout << "que push source vertex & source blossom elements" << endl;
	//cout << "order:";
	/*
	for (ll i = 0; i < order.size(); i++)
	{
		cout << order[i] << " ";
	}
	cout << endl;
	*/
	/////

	//-----step2-----
	while (!que.empty())
	{
		ll v = que.front();
		que.pop();

		//cout << "order:";
		/*
		for (ll i = 0; i < order.size(); i++)
		{
			cout << order[i] << " ";
		}
		cout << endl;
		*/

		ll v_idx;

		while (1)
		{
			ll u;
			ll mn_order = INF;
			ll mn_order_idx;

			/////
			/*
			cout << "que:";
			queue<ll> que2 = que;
			while (!que2.empty())
			{
				ll I = que2.front();
				que2.pop();
				cout << I << " ";
			}
			cout << endl;
			*/

			if (inBase[v])
			{
				for (ll i = 0; i < Ba.size(); i++)
				{
					if (Ba[i] == v)
					{
						v_idx = i;
						break;
					}
				}
				ll j = 0;
				for (ll i = 0; i < NBa.size(); i++)
				{
					if (j >= NBa_sub.size())
					{
						break;
					}
					if (NBa[i] == NBa_sub[j])
					{
						if (!C.X[v_idx][i].isZero())
						{
							if (p[NBa[i]] - p[v] == Q.X[v][NBa[i]])
							{
								if (order[NBa[i]] < INF && !judge_K(NBa[i], v, K))
								{
									//cout << K[NBa[i]] << " " << K[v] << endl;
									if (mn_order > order[NBa[i]])
									{
										u = NBa[i];
										mn_order = order[NBa[i]];
										mn_order_idx = i;
										////cout << "u:" << u << endl;
									}
								}
							}
						}
						j++;
					}
				}
			}
			else
			{
				for (ll i = 0; i < NBa.size(); i++)
				{
					if (NBa[i] == v)
					{
						v_idx = i;
						break;
					}
				}
				ll j = 0;
				for (ll i = 0; i < Ba.size(); i++)
				{
					if (j >= Ba_sub.size())
					{
						break;
					}
					if (Ba[i] == Ba_sub[j])
					{
						if (!C.X[i][v_idx].isZero())
						{
							if (p[v] - p[Ba[i]] == Q.X[v][Ba[i]])
							{
								if (order[Ba[i]] < INF && !judge_K(Ba[i], v, K))
								{
									//cout << K[Ba[i]] << " " << K[v] << endl;
									if (mn_order > order[Ba[i]])
									{
										u = Ba[i];
										mn_order = order[Ba[i]];
										mn_order_idx = i;
										//cout << "u:" << u << endl;
									}
								}
							}
						}
						j++;
					}
				}
			}
			if (mn_order == INF)
			{
				break;
			}

			//cout << v << " is adjacent to the labeled vertex:" << u << endl;
			//cout << "K[v]:" << K[v] << " K[u]:" << K[u] << endl;
			/*
			if (judge_K(u, v, K))
			{
				cout << "true" << endl;
			}
			else
			{
				cout << "false" << endl;
			}
			*/
			/*cout<<"P:"<<endl;
			for(ll i=0;i<P.size();i++)
			{
				if(!P[i].empty()){
					cout<<"P["<<i<<"]:";
					for(ll j=0;j<P[i].size();j++)
					{
						cout<<P[i][j]<<" ";
					}
					cout<<endl;
				}
			} */

			ll vs = v;
			while (P[vs][0] != vs)
			{
				vs = P[vs][0];
			}
			ll us = u;
			while (P[us][0] != us)
			{
				us = P[us][0];
			}
			if (vs != us && !isLine(vs, us))
			{
				//cout << endl;
				//cout << "-----step3-1------" << endl;
				path = Path(v, P);
				reverse(path.begin(), path.end());

				vector<ll> path_u = Path(u, P);

				/////
				/*
				cout << "P:";
				for (ll i = 0; i < P.size(); i++)
				{
					if (!P[i].empty())
					{
						for (ll j = 0; j < P[i].size(); j++)
						{
							cout << P[i][j] << " ";
						}
						cout << endl;
					}
				}
				cout << "RH";
				vector<node *> allblo = all_blossoms(Bl);
				for (ll i = 0; i < allblo.size(); i++)
				{
					node *Hi = allblo[i];
					cout << endl;
					if (Hi->normal)
					{
						cout << "normal" << endl;
					}
					else
					{
						cout << "source" << endl;
					}
					cout << "Hi key:" << Hi->key << endl;

					for (ll j = 0; j < Hi->ordering.size(); j++)
					{
						ll k = Hi->ordering[j];
						cout << k << " ";
						for (ll cnt = 0; cnt < Hi->routing[k].size(); cnt++)
						{
							cout << Hi->routing[k][cnt] << " ";
						}
						cout << endl;
					}
				}
				cout << "path(" << v << "):";
				for (ll cnt = 0; cnt < path.size(); cnt++)
				{
					cout << path[cnt] << " ";
				}
				cout << endl;
				cout << "path(" << u << "):";
				for (ll cnt = 0; cnt < path_u.size(); cnt++)
				{
					cout << path_u[cnt] << " ";
				}
				cout << endl;
				*/
				/////

				path.reserve(path.size() + path_u.size());
				path.insert(path.end(), path_u.begin(), path_u.end());

				return path;
			}
			else
			{
				//cout << endl;
				//cout << "-----step3-2------" << endl;
				Blossom(u, v, C, Q, P, Bl, Bl2, K, order, Ba, NBa, Ba_sub, NBa_sub, inBase, p, rho, que, N);
			}
		}

		vector<ll> unlabeled;
		if (inBase[v])
		{
			ll j = 0;
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (j >= NBa_sub.size())
				{
					break;
				}
				if (NBa[i] == NBa_sub[j])
				{
					if (!C.X[v_idx][i].isZero())
					{
						if (p[NBa[i]] - p[v] == Q.X[v][NBa[i]])
						{
							if (order[NBa[i]] == INF && rho[NBa[i]] == -1)
							{
								//cout << NBa[i] << " " << rho[NBa[i]] << endl;
								unlabeled.push_back(NBa[i]);
							}
						}
					}
					j++;
				}
			}
		}
		else
		{
			ll j = 0;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (j >= Ba_sub.size())
				{
					break;
				}
				if (Ba[i] == Ba_sub[j])
				{
					if (!C.X[i][v_idx].isZero())
					{
						if (p[v] - p[Ba[i]] == Q.X[v][Ba[i]])
						{
							if (order[Ba[i]] == INF && rho[Ba[i]] == -1)
							{

								//cout << Ba[i] << " " << rho[Ba[i]] << endl;
								unlabeled.push_back(Ba[i]);
							}
						}
					}
					j++;
				}
			}
		}
		/////
		/*
		cout << "unlabeled:";
		for (ll i = 0; i < unlabeled.size(); i++)
		{
			cout << unlabeled[i] << " ";
		}
		cout << endl;
		*/
		/////
		for (ll i = 0; i < unlabeled.size(); i++)
		{

			ll u = unlabeled[i];
			//cout << v<<" is adjacent to the unlabeled vertex " << u << endl;
			//cout << rho[u] << endl;
			if (order[u] < INF || rho[u] != -1)
			{
				continue;
			}

			if (isSingle(u, K))
			{
				//cout << endl;
				//cout << "-----step4-1------" << endl;

				P[mate(u)].push_back(v);
				P[mate(u)].push_back(u);
				order[mate(u)] = ordering_number;
				ordering_number++;
				rho[u] = v;
				que.push(mate(u));

				if (existsEdge(mate(u), v, C, Q, Ba, NBa, inBase, p))
				{
					Blossom(mate(u), v, C, Q, P, Bl, Bl2, K, order, Ba, NBa, Ba_sub, NBa_sub, inBase, p, rho, que, N);
				}
			}
			else
			{
				node *Hi = tree_search(Bl, K[u]).first;
				if (Hi->normal)
				{
					ll bi = Hi->bud;

					if (existsEdge(v, bi, C, Q, Ba, NBa, inBase, p))
					{
						//cout << endl;
						//cout << "-----step4-2------" << endl;
						if (u != bi)
						{
							//cout << "-----Graft start-----" << endl;
							Graft(v, Hi, C, Q, P, Bl, Bl2, K, Ba, NBa, Ba_sub, NBa_sub, inBase, order, p, rho, que);
							//cout << "------Graft end-------" << endl;
						}
						else
						{
							P[Hi->tip] = {v, bi};
							for (ll j = 0; j < order.size(); j++)
							{
								if (order[j] < INF && order[j] > order[v])
								{
									order[j]++;
								}
							}
							order[Hi->tip] = order[v] + 1;
							ordering_number++;

							rho[bi] = v;
							Hi->label = 1;
							for (ll j = 0; j < Hi->ordering.size(); j++)
							{
								if (Hi->ordering[j] != Hi->tip)
								{
									P[Hi->ordering[j]] = Hi->routing[Hi->ordering[j]];
									que.push(Hi->ordering[j]);
									order[Hi->ordering[j]] = ordering_number;
									ordering_number++;
								}
							}
						}
					}
					else
					{
						//cout << endl;
						//cout << "-----step4-3------" << endl;
						//cout << "K[u]:" << K[u] << " " << Hi->key << endl;
						ll y;
						for (ll j = 0; j < Hi->ordering.size(); j++)
						{
							if (existsEdge(v, Hi->ordering[j], C, Q, Ba, NBa, inBase, p))
							{
								y = Hi->ordering[j];
								break;
							}
						}

						Hi->label = -1;
						P[bi].push_back(v);
						order[bi] = ordering_number;
						ordering_number++;
						ll k = y;
						P[bi].push_back(y);
						while (Hi->routing[k][0] != k)
						{
							for (ll j = Hi->routing[k].size() - 1; j >= 0; j--)
							{
								P[bi].push_back(Hi->routing[k][j]);
							}
							k = Hi->routing[k][0];
						}
						que.push(bi);

						for (ll j = 0; j < Hi->ordering.size(); j++)
						{
							if (order[Hi->ordering[j]] == INF)
							{
								rho[Hi->ordering[j]] = v;
							}
						}
					}
				}
			}
		}
	}
	////cout << "-----------Search end------------------" << endl;
	return {};
}

bool DualUpdate(Matrix &C, Tree &Bl, vector<Field> &p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &order, vector<ll> &rho, vector<ll> &K, vector<bool> &inBase, Matrix &Q)
{
	if (!overflow)
	{
		DualUpdate_num++;
	}
	//cout << "DualUpdate start" << endl;
	/*
	for (ll i = 0; i < Ba.size(); i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			ll v = Ba[i];
			ll u = NBa[j];

			ll v_idx, u_idx;
			for (ll k = 0; k < Ba.size(); k++)
			{
				if (Ba[k] == v)
				{
					v_idx = k;
					break;
				}
			}
			for (ll k = 0; k < NBa.size(); k++)
			{
				if (NBa[k] == u)
				{
					u_idx = k;
					break;
				}
			}
			if (!C.X[v_idx][u_idx].isZero())
			{
				if (p[u] - p[v] == Q.X[u][v])
				{
					cout << "tight:" << v << " " << u << endl;
				}
				else
				{
					cout << "no tight:" << v << " " << u << endl;
				}
			}
		}
	}
	*/

	vector<ll> R_plus;
	vector<ll> R_minus;
	vector<ll> Z_plus;
	vector<ll> Z_minus;

	ll N = rho.size();
	vector<bool> flag;
	flag.resize(N);

	for (int i = 0; i < N; i++)
	{
		flag[i] = false;
	}

	for (ll i = 0; i < Ba.size(); i++)
	{
		ll k = Ba[i];
		if (order[k] < INF)
		{
			if (isSingle(k, K))
			{
				R_plus.push_back(k);
				flag[k] = true;
			}
			else
			{
				node *H = tree_search(Bl, K[k]).first;
				if (H->normal)
				{
					if (k == H->bud)
					{
						R_plus.push_back(k);
						flag[k] = true;
					}
				}
			}
		}
		else if (rho[k] != -1)
		{
			if (isSingle(k, K))
			{
				R_minus.push_back(k);
				flag[k] = true;
			}
			else
			{
				node *H = tree_search(Bl, K[k]).first;
				if (H->normal)
				{
					if (k == H->bud)
					{
						R_minus.push_back(k);
						flag[k] = true;
					}
				}
			}
		}
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		ll k = NBa[i];
		if (order[k] < INF)
		{
			if (isSingle(k, K))
			{
				R_plus.push_back(k);
				flag[k] = true;
			}
			else
			{
				node *H = tree_search(Bl, K[k]).first;
				if (H->normal)
				{
					if (k == H->bud)
					{
						R_plus.push_back(k);
						flag[k] = true;
					}
				}
			}
		}
		else if (rho[k] != -1)
		{
			if (isSingle(k, K))
			{
				R_minus.push_back(k);
				flag[k] = true;
			}
			else
			{
				node *H = tree_search(Bl, K[k]).first;
				if (H->normal)
				{
					if (k == H->bud)
					{
						R_minus.push_back(k);
						flag[k] = true;
					}
				}
			}
		}
	}

	vector<node *> MaximalBlossom = maximal_blossom(Bl);

	for (ll i = 0; i < MaximalBlossom.size(); i++)
	{
		node *Hi = MaximalBlossom[i];
		if (Hi->label == 1)
		{
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				Z_plus.push_back(Hi->ordering[j]);
				flag[Hi->ordering[j]] = true;
			}
		}
		if (Hi->label == -1)
		{
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				Z_minus.push_back(Hi->ordering[j]);
				flag[Hi->ordering[j]] = true;
			}
		}
	}

	//Yを定める
	vector<ll> Y;
	for (ll i = 0; i < Ba.size(); i++)
	{
		if (flag[Ba[i]] == false)
		{
			Y.push_back(Ba[i]);
		}
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		if (flag[NBa[i]] == false)
		{
			Y.push_back(NBa[i]);
		}
	}

	//define eps
	set<ll> RZ_plus;
	Field eps1, eps2, eps3, eps4;
	for (ll i = 0; i < R_plus.size(); i++)
	{
		RZ_plus.insert(R_plus[i]);
	}
	for (ll i = 0; i < Z_plus.size(); i++)
	{
		RZ_plus.insert(Z_plus[i]);
	}

	/////
	/*
	cout << "R+:";
	for (ll i = 0; i < R_plus.size(); i++)
	{
		cout << R_plus[i] << " ";
	}
	cout << endl;
	cout << "R-:";
	for (ll i = 0; i < R_minus.size(); i++)
	{
		cout << R_minus[i] << " ";
	}
	cout << endl;
	cout << "Z+:";
	for (ll i = 0; i < Z_plus.size(); i++)
	{
		cout << Z_plus[i] << " ";
	}
	cout << endl;
	cout << "Z-:";
	for (ll i = 0; i < Z_minus.size(); i++)
	{
		cout << Z_minus[i] << " ";
	}
	cout << endl;
	cout << "Y:";
	for (ll i = 0; i < Y.size(); i++)
	{
		cout << Y[i] << " ";
	}
	cout << endl;
	*/
	/////

	vector<ll> V_inverse;
	V_inverse.resize(order.size());

	for (ll i = 0; i < Ba.size(); i++)
	{
		V_inverse[Ba[i]] = i;
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		V_inverse[NBa[i]] = i;
	}
	//eps1
	int ope1 = 0;
	for (auto itr1 = RZ_plus.begin(); itr1 != RZ_plus.end(); itr1++)
	{
		for (auto itr2 = RZ_plus.begin(); itr2 != RZ_plus.end(); itr2++)
		{
			if (inBase[*itr1] == 1 && inBase[*itr2] == 0)
			{
				if (!C.X[V_inverse[*itr1]][V_inverse[*itr2]].isZero() && !judge_K(*itr1, *itr2, K))
				{
					if (ope1 == 0)
					{
						eps1 = (p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
						ope1 = 1;
					}
					else
					{
						eps1 = min(eps1, p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
					}
				}
			}
			else if (inBase[*itr1] == 0 && inBase[*itr2] == 1)
			{
				if (!C.X[V_inverse[*itr2]][V_inverse[*itr1]].isZero() && !judge_K(*itr1, *itr2, K))
				{
					if (ope1 == 0)
					{
						eps1 = p[*itr1] - p[*itr2] - Q.X[*itr1][*itr2];
						ope1 = 1;
					}
					else
					{
						eps1 = min(eps1, p[*itr1] - p[*itr2] - Q.X[*itr1][*itr2]);
					}
				}
			}
		}
	}
	if (ope1 != 0)
	{
		eps1 = eps1.divideByTwo();
	}

	//eps2
	int ope2 = 0;
	for (auto itr1 = RZ_plus.begin(); itr1 != RZ_plus.end(); itr1++)
	{
		for (auto itr2 = Y.begin(); itr2 != Y.end(); itr2++)
		{
			if (inBase[*itr1] == 1 && inBase[*itr2] == 0)
			{
				if (!C.X[V_inverse[*itr1]][V_inverse[*itr2]].isZero() && !judge_K(*itr1, *itr2, K))
				{
					if (ope2 == 0)
					{
						eps2 = (p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
						ope2 = 1;
					}
					else
					{
						eps2 = min(eps2, p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
					}
				}
			}
		}
	}

	//eps3
	int ope3 = 0;
	for (auto itr1 = Y.begin(); itr1 != Y.end(); itr1++)
	{
		for (auto itr2 = RZ_plus.begin(); itr2 != RZ_plus.end(); itr2++)
		{
			if (inBase[*itr1] == 1 && inBase[*itr2] == 0)
			{
				if (!C.X[V_inverse[*itr1]][V_inverse[*itr2]].isZero() && !judge_K(*itr1, *itr2, K))
				{
					if (ope3 == 0)
					{
						eps3 = (p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
						ope3 = 1;
					}
					else
					{
						eps3 = min(eps3, p[*itr2] - p[*itr1] - Q.X[*itr1][*itr2]);
					}
				}
			}
		}
	}

	//eps4
	int ope4 = 0;
	for (ll i = 0; i < MaximalBlossom.size(); i++)
	{
		node *Hi = MaximalBlossom[i];
		if (Hi->label == -1)
		{
			if (ope4 == 0)
			{
				eps4 = Hi->q;
			}
			else
			{
				eps4 = min(eps4, Hi->q);
			}
			ope4 = 1;
		}
	}

	/////
	/*
	if (ope1 == 1)
	{
		cout << "eps1:";
		eps1.output();
		cout << endl;
	}
	else
	{
		cout << "eps1 is INF" << endl;
	}
	if (ope2 == 1)
	{
		cout << "eps2:";
		eps2.output();
		cout << endl;
	}
	else
	{
		cout << "eps2 is INF" << endl;
	}
	if (ope3 == 1)
	{
		cout << "eps3:";
		eps3.output();
		cout << endl;
	}
	else
	{
		cout << "eps3 is INF" << endl;
	}
	if (ope4 == 1)
	{
		cout << "eps4:";
		eps4.output();
		cout << endl;
	}
	else
	{
		cout << "eps4 is INF" << endl;
	}
	*/
	/////

	Field eps;
	if (ope1 + ope2 + ope3 + ope4 == 0)
	{
		//cout << "DualUpdate end" << endl;
		return false;
	}
	else if (ope2 + ope3 + ope4 == 0)
	{
		eps = eps1;
	}
	else if (ope1 + ope3 + ope4 == 0)
	{
		eps = eps2;
	}
	else if (ope1 + ope2 + ope4 == 0)
	{
		eps = eps3;
	}
	else if (ope1 + ope2 + ope3 == 0)
	{
		eps = eps4;
	}
	else if (ope3 + ope4 == 0)
	{
		eps = min(eps1, eps2);
	}
	else if (ope2 + ope4 == 0)
	{
		eps = min(eps1, eps3);
	}
	else if (ope2 + ope3 == 0)
	{
		eps = min(eps1, eps4);
	}
	else if (ope1 + ope4 == 0)
	{
		eps = min(eps2, eps3);
	}
	else if (ope1 + ope3 == 0)
	{
		eps = min(eps2, eps4);
	}
	else if (ope1 + ope2 == 0)
	{
		eps = min(eps3, eps4);
	}
	else if (ope1 == 0)
	{
		eps = min(eps2, min(eps3, eps4));
	}
	else if (ope2 == 0)
	{
		eps = min(eps1, min(eps3, eps4));
	}
	else if (ope3 == 0)
	{
		eps = min(eps1, min(eps2, eps4));
	}
	else if (ope4 == 0)
	{
		eps = min(eps1, min(eps2, eps3));
	}
	else
	{
		eps = min(min(eps1, eps2), min(eps3, eps4));
	}

	/*
	cout << "eps:";
	eps.output();
	cout << endl;
	*/
	Field q;
	q.setZero();
	if (eps < q)
	{
		exit(EXIT_FAILURE);
	}

	//SQ.output_matrix();

	//update dual variables
	//cout<<"x"<<endl;
	for (ll i = 0; i < R_plus.size(); i++)
	{
		ll v = R_plus[i];
		if (inBase[v] == 1)
		{
			p[v] += eps;
		}
		else
		{
			p[v] -= eps;
		}
	}
	for (ll i = 0; i < R_minus.size(); i++)
	{
		ll v = R_minus[i];
		if (inBase[v] == 1)
		{
			p[v] -= eps;
		}
		else
		{
			p[v] += eps;
		}
	}
	//cout<<"y"<<endl;

	for (ll i = 0; i < MaximalBlossom.size(); i++)
	{
		node *Hi = MaximalBlossom[i];

		if (Hi->label == 1)
		{
			Hi->q += eps;

			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				ll v = Hi->ordering[j];
				for (ll k = 0; k < Ba.size(); k++)
				{
					Q.X[v][Ba[k]] += eps;
					Q.X[Ba[k]][v] += eps;
					//cout<<k<<endl;
				}
				for (ll k = 0; k < NBa.size(); k++)
				{
					Q.X[v][NBa[k]] += eps;
					Q.X[NBa[k]][v] += eps;
					//cout<<k<<endl;
				}
			}
			//cout << "a" << endl;
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				for (ll k = 0; k < Hi->ordering.size(); k++)
				{
					ll u = Hi->ordering[j];
					ll v = Hi->ordering[k];
					if (u != v)
					{
						Q.X[u][v] = Q.X[u][v] - eps;
						Q.X[v][u] = Q.X[v][u] - eps;
					}
					else
					{ //Q_uuも二回ずつ加算されている
						Q.X[u][u] -= eps;
					}
				}
			}
			//cout << "b" << endl;
		}
		else if (Hi->label == -1)
		{
			Hi->q -= eps;
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				ll v = Hi->ordering[j];
				for (ll k = 0; k < Ba.size(); k++)
				{
					Q.X[v][Ba[k]] -= eps;
					Q.X[Ba[k]][v] -= eps;
				}
				for (ll k = 0; k < NBa.size(); k++)
				{
					Q.X[v][NBa[k]] -= eps;
					Q.X[NBa[k]][v] -= eps;
				}
			}
			////cout << "c" << endl;
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				for (ll k = 0; k < Hi->ordering.size(); k++)
				{
					ll u = Hi->ordering[j];
					ll v = Hi->ordering[k];
					if (u != v)
					{
						Q.X[u][v] += eps;
						Q.X[v][u] += eps;
					}
					else
					{
						Q.X[u][u] += eps;
					}
				}
			}
			//cout << "d" << endl;
		}
	}

	/*
	cout << "p:";
	for (ll i = 0; i < p.size(); i++)
	{
		p[i].output();
		cout << " ";
	}
	cout << endl;
	vector<node *> MaxBlo = maximal_blossom(Bl);
	for (ll i = 0; i < MaxBlo.size(); i++)
	{
		cout << "q(H" << MaxBlo[i]->key << "):";
		MaxBlo[i]->q.output();
		cout << endl;
	}

	cout << "Q:";
	Q.output_matrix();
	cout << endl;
	cout << endl;
	Q.X[9][20].output();
	cout << "DualUpdate end" << endl;
	*/

	return true;
}

void Augment(vector<ll> &path, Tree &Bl, Matrix &C, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<vector<ll>> &P, vector<ll> &K, vector<ll> &order, vector<ll> &rho, vector<bool> &inBase, Matrix &Q, vector<Field> &p)
{
	if (!overflow)
	{
		Augment_num++;
	}
	//cout << "-------Augment start-----------" << endl;

	/*
	vector<node *> all_blossom = all_blossoms(Bl);
	for (ll i = 0; i < all_blossom.size(); i++)
	{
		node *H_ = all_blossom[i];
	//	cout << "key:" << H_->key << " "
	//		 << "elements:";
	/*
		for (ll j = 0; j < H_->ordering.size(); j++)
		{
			cout << H_->ordering[j] << " ";
		}
		cout << endl;
		cout << "q:";
		H_->q.output();
		cout << endl;
	}
	*/
	//-----step0-----
	vector<node *> L = all_blossoms(Bl);
	//cout << "C:";
	//C.output_matrix();
	/*
	cout << endl;
	cout << "path:";
	for (ll i = 0; i < path.size(); i++)
	{
		cout << path[i] << " ";
	}
	cout << endl;
	*/

	set<ll> Lp;
	vector<node *> NLp;
	vector<pair<node *, pair<ll, ll>>> Lp_positive;
	set<node *> Lp_zero;

	//Λpを求める
	vector<vector<ll>> S;
	S.resize(rho.size());
	for (ll i = 0; i < L.size(); i++)
	{
		node *Hi = L[i];
		for (ll j = 0; j < Hi->ordering.size(); j++)
		{
			ll v = Hi->ordering[j];
			S[v].push_back(Hi->key);
		}
	}
	for (ll i = 0; i < S.size(); i++)
	{
		sort(S[i].begin(), S[i].end());
	}

	for (ll i = 0; i < path.size() / 2; i++)
	{
		ll u = path[2 * i];
		ll v = path[2 * i + 1];

		vector<ll> sym_dif;
		sym_dif.resize(S[u].size() + S[v].size());
		auto itr = set_symmetric_difference(S[u].begin(), S[u].end(), S[v].begin(), S[v].end(), sym_dif.begin());
		sym_dif.resize(itr - sym_dif.begin());
		for (ll j = 0; j < sym_dif.size(); j++)
		{
			node *H = tree_search(Bl, sym_dif[j]).first;
			Lp.insert(H->key);
			if (H->q.isZero())
			{
				Lp_zero.insert(H);
			}
			else
			{
				Lp_positive.push_back(make_pair(H, make_pair(u, v)));
			}
		}
	}

	//cout << "all blossom:";
	for (ll i = 0; i < L.size(); i++)
	{
		node *Hi = L[i];
		//cout << L[i]->key << " ";

		if (Lp.find(Hi->key) == Lp.end())
		{
			NLp.push_back(Hi);
		}
	}
	//cout << endl;
	//cout << "Lp_positive:" << endl;
	/*
	for (ll i = 0; i < Lp_positive.size(); i++)
	{
		cout << Lp_positive[i].first->key << " " << Lp_positive[i].second.first << " " << Lp_positive[i].second.second << endl;
	}
	cout << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!Lp_zero";
	*/
	/*
	for (auto itr = Lp_zero.begin(); itr != Lp_zero.end(); itr++)
	{
		cout << (*itr)->key << " ";
	}
	cout << endl;
	*/

	//!!!!!変えた!!!!! isMaximal(Bl,H)

	for (ll i = 0; i < NLp.size(); i++)
	{
		node *H = NLp[i];
		if (H->q.isZero() && isMaximal(Bl, H))
		{
			Expand(C, H, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Bl);
		}
	}

	//-----step1-----

	vector<pr> M;
	//MにPのマッチング(Cのいんでっくすのじょうたいで)突っ込む
	for (ll i = 0; i < path.size() / 2; i++)
	{
		pr p;
		ll u = path[2 * i];
		ll v = path[2 * i + 1];
		if (inBase[u] == 1 && inBase[v] == 0)
		{
			ll idx_u, idx_v;
			for (ll j = 0; j < Ba.size(); j++)
			{
				if (Ba[j] == u)
				{
					idx_u = j;
				}
			}
			for (ll j = 0; j < NBa.size(); j++)
			{
				if (NBa[j] == v)
				{
					idx_v = j;
				}
			}
			p = pr(idx_u, idx_v);
		}
		if (inBase[u] == 0 && inBase[v] == 1)
		{
			ll idx_u, idx_v;
			for (ll j = 0; j < NBa.size(); j++)
			{
				if (NBa[j] == u)
				{
					idx_u = j;
				}
			}
			for (ll j = 0; j < Ba.size(); j++)
			{
				if (Ba[j] == v)
				{
					idx_v = j;
				}
			}
			p = pr(idx_v, idx_u);
		}
		M.push_back(p);
	}

	vector<pr> new_dummy_lines;
	for (ll I = 0; I < Lp_positive.size(); I++)
	{
		node *H = Lp_positive[I].first;
		ll u = Lp_positive[I].second.first;
		ll v = Lp_positive[I].second.second;

		ll xi;
		ll yi;
		for (ll j = 0; j < H->ordering.size(); j++)
		{
			if (u == H->ordering[j])
			{
				xi = u;
				yi = v;
				break;
			}
			if (v == H->ordering[j])
			{
				xi = v;
				yi = u;
				break;
			}
		}

		//~~~~~xiとyiを求める
		ll ti2 = new_vertex(K, P, C, Q, order, rho, inBase, p);
		ll bi2 = new_vertex(K, P, C, Q, order, rho, inBase, p);

		//cout << bi2 << " " << ti2 << " " << xi << " " << yi << endl;

		new_dummy_lines.push_back(pr(ti2, bi2));

		vector<node *> H_ancestors = ancestors(Bl, H);
		H->ordering.push_back(ti2);
		//H->ordering.push_back(bi2); //あとで消す
		set<ll> H_elements;
		for (ll i = 0; i < H->ordering.size() - 1; i++)
		{
			H_elements.insert(H->ordering[i]);
		}

		Field q_ancestors;
		q_ancestors.setZero();
		for (ll i = 0; i < H_ancestors.size(); i++)
		{
			node *Hj = H_ancestors[i];
			Hj->ordering.push_back(bi2);
			Hj->ordering.push_back(ti2);
			q_ancestors += Hj->q;
		}

		if (H->normal)
		{
			for (ll i = 0; i < Ba.size(); i++)
			{
				Q.X[Ba[i]][ti2] = Q.X[Ba[i]][H->tip];
				Q.X[ti2][Ba[i]] = Q.X[H->tip][Ba[i]];

				Q.X[Ba[i]][bi2] = Q.X[Ba[i]][H->bud];
				Q.X[bi2][Ba[i]] = Q.X[H->bud][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][ti2] = Q.X[NBa[i]][H->tip];
				Q.X[ti2][NBa[i]] = Q.X[H->tip][NBa[i]];

				Q.X[NBa[i]][bi2] = Q.X[NBa[i]][H->bud];
				Q.X[bi2][NBa[i]] = Q.X[H->bud][NBa[i]];
			}
			Q.X[ti2][H->tip].setZero();
			Q.X[H->tip][ti2].setZero();
			Q.X[bi2][H->bud].setZero();
			Q.X[H->bud][bi2].setZero();
			Q.X[ti2][ti2] = Q.X[H->tip][H->tip];
			Q.X[bi2][bi2] = Q.X[H->bud][H->bud];
			Q.X[ti2][bi2] = H->q;
			Q.X[bi2][ti2] = H->q;
		}
		else
		{
			//cout << "!!!!!!!!!!!!!!!!!!!!!!" << endl;
			//cout << "H->q:";
			//H->q.output();
			//cout << endl;
			//cout << "q_ancestors:";
			//q_ancestors.output();
			//cout << endl;
			Q.X[ti2][ti2] = q_ancestors + H->q;
			Q.X[bi2][bi2] = q_ancestors;
			for (ll i = 0; i < Ba.size(); i++)
			{
				ll k = Ba[i];
				Q.X[ti2][k] += Q.X[k][k] + Q.X[ti2][ti2];
				Q.X[k][ti2] += Q.X[k][k] + Q.X[ti2][ti2];
				Q.X[bi2][k] += Q.X[k][k] + Q.X[bi2][bi2];
				Q.X[k][bi2] += Q.X[k][k] + Q.X[bi2][bi2];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				ll k = NBa[i];
				Q.X[ti2][k] += Q.X[k][k] + Q.X[ti2][ti2];
				Q.X[k][ti2] += Q.X[k][k] + Q.X[ti2][ti2];
				Q.X[bi2][k] += Q.X[k][k] + Q.X[bi2][bi2];
				Q.X[k][bi2] += Q.X[k][k] + Q.X[bi2][bi2];
			}

			for (ll i = 0; i < H->ordering.size(); i++)
			{
				ll k = H->ordering[i];
				if (k != ti2)
				{
					Q.X[ti2][k] -= H->q.Double();
					Q.X[k][ti2] -= H->q.Double();
				}
			}

			for (ll i = 0; i < H_ancestors.size(); i++)
			{
				node *Hi = H_ancestors[i];
				for (ll j = 0; j < Hi->ordering.size(); j++)
				{
					ll k = Hi->ordering[j];
					if (k != ti2 && k != bi2)
					{
						Q.X[ti2][k] -= Hi->q.Double();
						Q.X[k][ti2] -= Hi->q.Double();
						Q.X[bi2][k] -= Hi->q.Double();
						Q.X[k][bi2] -= Hi->q.Double();
					}
				}
			}
			Q.X[ti2][bi2] = H->q;
			Q.X[bi2][ti2] = H->q;
		}

		if (inBase[xi] == 1 && inBase[yi] == 0)
		{
			Ba.push_back(bi2);
			NBa.push_back(ti2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[bi2] = true;
			inBase[ti2] = false;

			//Blossom step2 で　b<-bi t<-ti r<-xi g<-yi (Hに入ってるか入ってないかなど少し違うので注意)
			ll bi2_idx = Ba.size() - 1;
			ll ti2_idx = NBa.size() - 1;
			ll xi_idx;
			ll yi_idx;

			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == xi)
				{
					xi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == yi)
				{
					yi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(NBa[i]) == H_elements.end())
				{
					C.X[bi2_idx][i] = C.X[xi_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (H_elements.find(Ba[i]) != H_elements.end())
				{
					C.X[i][ti2_idx] = C.X[i][yi_idx];
				}
			}
			C.X[bi2_idx][ti2_idx].setZero();

			p[bi2] = p[yi] - Q.X[bi2][yi];
			p[ti2] = p[xi] + Q.X[xi][ti2];

			M.push_back(pr(bi2_idx, ti2_idx));
		}
		else if (inBase[xi] == 0 && inBase[yi] == 1)
		{
			Ba.push_back(ti2);
			NBa.push_back(bi2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[bi2] = false;
			inBase[ti2] = true;

			ll ti2_idx = Ba.size() - 1;
			ll bi2_idx = NBa.size() - 1;
			ll xi_idx;
			ll yi_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == yi)
				{
					yi_idx = i;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == xi)
				{
					xi_idx = i;
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				//if NBa[i] in H\B*
				if (H_elements.find(Ba[i]) == H_elements.end())
				{
					C.X[i][bi2_idx] = C.X[i][xi_idx];
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				if (H_elements.find(NBa[i]) != H_elements.end())
				{
					C.X[ti2_idx][i] = C.X[yi_idx][i];
				}
			}
			C.X[ti2_idx][bi2_idx].setZero();

			p[bi2] = p[yi] + Q.X[bi2][yi];
			p[ti2] = p[xi] - Q.X[xi][ti2];

			M.push_back(pr(ti2_idx, bi2_idx));
		}
	}

	//-----step2-----
	Pivoting_around_M(C, M, Ba, NBa, Ba_sub, NBa_sub, inBase);
	/////
	/*
	if(!Pivoting_around_M(C, M, Ba, NBa, Ba_sub, NBa_sub, inBase)){
		cout<<"path/K"<<endl;
		for(ll i=0;i<path.size();i++){
			cout<<"v:"<<path[i]<<" K[v]:"<<K[path[i]]<<" p[v]:";p[path[i]].output();cout<<endl;;
		}
		cout<<"Lp+:"<<endl;
		for(ll i=0;i<Lp_positive.size();i++){
			cout<<Lp_positive[i].first->key<<" "<<Lp_positive[i].second.first<<" "<<Lp_positive[i].second.second;
		}
		cout<<"new dummy line"<<endl;
		for(ll i=0;i<new_dummy_lines.size();i++){
			cout<<new_dummy_lines[i].first<<" "<<new_dummy_lines[i].second<<endl;
		}

		for(ll i=0;i<M.size();i++){
			ll v=Ba[M[i].first];
			for(ll j=0;j<M.size();j++)
			{
				ll u=NBa[M[j].second];
				
				if(!C.X[M[i].first][M[j].second].isZero())
				{
					if(p[u]-p[v]==Q.X[u][v]){
						cout<<"tight:"<<v<<" "<<u<<endl;
					}
					else{
						cout<<"not tight:"<<v<<" "<<u<<endl;
					}
				}
			}
		}

		cout<<"M:"<<endl;
		for(ll i=0;i<M.size();i++){
			cout<<M[i].first<<" "<<M[i].second<<" "<<Ba[M[i].first]<<" "<<NBa[M[i].second]<<endl;
		}
		
		exit(EXIT_FAILURE);
		
	}
	*/

	//-----step3-----

	//cout << "step3" << endl;
	for (auto itr_z = Lp_zero.begin(); itr_z != Lp_zero.end(); itr_z++)
	{
		node *Hi = *itr_z;
		//cout << "Hi->key:" << Hi->key << endl;
		/*
		if (Hi->normal)
		{
			cout << "normal" << endl;
		}
		else
		{
			cout << "source" << endl;
		}
		*/

		if (!Hi->normal)
		{
			tree_delete(Bl, Hi->key);
		}
		else
		{
			ll ti = Hi->tip;
			ll bi = Hi->bud;

			vector<node *> H_ancestors = ancestors(Bl, Hi);
			for (ll i = 0; i < H_ancestors.size(); i++)
			{
				//cout << "sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
				node *Hj = H_ancestors[i];
				ll idx;
				for (ll j = 0; j < Hj->ordering.size(); j++)
				{
					ll k = Hj->ordering[j];
					if (k == ti)
					{
						idx = j;
						break;
					}
				}
				auto itr = Hj->ordering.erase(Hj->ordering.begin() + idx);

				for (ll j = 0; j < Hj->ordering.size(); j++)
				{
					ll k = Hj->ordering[j];
					if (k == bi)
					{
						idx = j;
						break;
					}
				}
				itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
			}
			tree_delete(Bl, Hi->key);
			remove_vertex(ti, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
			remove_vertex(bi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
		}
		//if (Hi->normal)
		//{
		//Expand(C, Hi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Bl);
		//}
	}
	//cout << "C:";
	//C.output_matrix();
	//cout << endl;

	//rename bi',ti' as bud tip of Hi

	for (ll I = 0; I < Lp_positive.size(); I++)
	{
		node *Hi = Lp_positive[I].first;
		//cout << Hi->key << endl;
		if (Hi->normal)
		{
			ll former_tip = Hi->tip;
			ll former_bud = Hi->bud;
			remove_vertex(Hi->bud, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
			//cout << "eee" << endl;
			remove_vertex(Hi->tip, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);

			//cout << "aaa" << endl;
			Hi->tip = new_dummy_lines[I].first;
			Hi->bud = new_dummy_lines[I].second;

			ll idx;
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				if (Hi->ordering[j] == former_tip)
				{
					idx = j;
					break;
				}
			}
			Hi->ordering.erase(Hi->ordering.begin() + idx);
			//Hi->tip = Hi->ordering[Hi->ordering.size() - 1];
			//Hi->bud = Hi->ordering[Hi->ordering.size() - 1];
			//Hi->ordering.erase(Hi->ordering.end() - 1);
			//Hi->ordering.erase(Hi->ordering.begin());

			//cout << "bbb" << endl;
			//cout << "new tip is " << Hi->tip << " /new bud is " << Hi->bud << endl;

			vector<node *> Hi_ancestors = ancestors(Bl, Hi);
			for (ll j = 0; j < Hi_ancestors.size(); j++)
			{
				node *Hj = Hi_ancestors[j];
				ll idx;
				for (ll k = 0; k < Hj->ordering.size(); k++)
				{
					if (Hj->ordering[k] == former_tip)
					{
						idx = k;
						break;
					}
				}

				auto itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
				for (ll k = 0; k < Hj->ordering.size(); k++)
				{
					if (Hj->ordering[k] == former_bud)
					{
						idx = k;
						break;
					}
				}
				itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
			}
		}
		else
		{
			/* Hi->tip = Hi->ordering[Hi->ordering.size() - 2];
			Hi->bud = Hi->ordering[Hi->ordering.size() - 1]; */
			Hi->tip = new_dummy_lines[I].first;
			Hi->bud = new_dummy_lines[I].second;
			//cout << "Key:" << Hi->key << endl;
			//cout << "new tip is " << Hi->tip << " new bud is:" << Hi->bud << endl;
			//Hi->ordering.erase(Hi->ordering.end() - 1);
			Hi->normal = true;
		}
	}

	//-----step4-----
	//cout << "-----step4-----" << endl;
	vector<node *> Blossoms = all_blossoms(Bl);
	vector<node *> NormalBlossoms;
	for (ll i = 0; i < Blossoms.size(); i++)
	{
		node *Hi = Blossoms[i];
		if (Hi->normal)
		{
			NormalBlossoms.push_back(Hi);
		}
	}

	vector<node *> Blossoms_in_MaximalNormalBlossom = blossoms_in_MaximalNormalBlossom(Bl);

	//cout << "!!!!!!Blossoms_in_MaximalNormalBlossom" << endl;
	/*
	for (ll i = 0; i < Blossoms_in_MaximalNormalBlossom.size(); i++)
	{
		node *Hi = Blossoms_in_MaximalNormalBlossom[i];
		//cout << "key:" << Hi->key << " "
		//	 << "elements:";
		for (ll j = 0; j < Hi->ordering.size(); j++)
		{
			cout << Hi->ordering[j] << " ";
		}
		cout << endl;
	}
	*/

	for (ll I = Blossoms_in_MaximalNormalBlossom.size() - 1; I >= 0; I--)
	{
		//node *Hi = NormalBlossoms[I];
		node *Hi = Blossoms_in_MaximalNormalBlossom[I];
		ll ti2 = new_vertex(K, P, C, Q, order, rho, inBase, p);
		ll bi2 = new_vertex(K, P, C, Q, order, rho, inBase, p);

		vector<node *> Hi_ancestors = ancestors(Bl, Hi);

		Hi->ordering.push_back(ti2);
		Hi->ordering.push_back(bi2);
		for (ll j = 0; j < Hi_ancestors.size(); j++)
		{
			node *Hj = Hi_ancestors[j];
			Hj->ordering.push_back(ti2);
			Hj->ordering.push_back(bi2);
		}
		set<ll> Hi_elements;
		for (ll j = 0; j < Hi->ordering.size() - 1; j++)
		{
			Hi_elements.insert(Hi->ordering[j]);
		}
		ll ti = Hi->tip;
		ll bi = Hi->bud;

		vector<pr> M2;
		if (inBase[bi] && !inBase[ti])
		{
			Ba.push_back(ti2);
			NBa.push_back(bi2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[ti2] = true;
			inBase[bi2] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				Q.X[Ba[i]][ti2] = Q.X[Ba[i]][ti];
				Q.X[ti2][Ba[i]] = Q.X[ti][Ba[i]];

				Q.X[Ba[i]][bi2] = Q.X[Ba[i]][bi];
				Q.X[bi2][Ba[i]] = Q.X[bi][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][ti2] = Q.X[NBa[i]][ti];
				Q.X[ti2][NBa[i]] = Q.X[ti][NBa[i]];

				Q.X[NBa[i]][bi2] = Q.X[NBa[i]][bi];
				Q.X[bi2][NBa[i]] = Q.X[bi][NBa[i]];
			}
			Q.X[ti2][ti].setZero();
			Q.X[ti][ti2].setZero();
			Q.X[bi2][bi].setZero();
			Q.X[bi][bi2].setZero();
			Q.X[ti2][ti2] = Q.X[ti][ti];
			Q.X[bi2][bi2] = Q.X[bi][bi];
			Q.X[ti2][bi2] = Hi->q;
			Q.X[bi2][ti2] = Hi->q;

			ll ti2_idx = Ba.size() - 1;
			ll bi2_idx = NBa.size() - 1;
			ll bi_idx;
			ll ti_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == bi)
				{
					bi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == ti)
				{
					ti_idx = i;
					break;
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in Hi\B*
				if (Hi_elements.find(NBa[i]) != Hi_elements.end())
				{
					C.X[ti2_idx][i] = C.X[bi_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (Hi_elements.find(Ba[i]) == Hi_elements.end())
				{
					C.X[i][bi2_idx] = C.X[i][ti_idx];
				}
			}

			p[bi2] = p[ti] - Q.X[bi2][ti];
			p[ti2] = p[bi] + Q.X[bi2][ti];

			M2.push_back(pr(bi_idx, ti_idx));
			M2.push_back(pr(ti2_idx, bi2_idx));
		}
		else if (inBase[bi] == 0 && inBase[ti] == 1)
		{
			Ba.push_back(bi2);
			NBa.push_back(ti2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[bi2] = true;
			inBase[ti2] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				Q.X[Ba[i]][ti2] = Q.X[Ba[i]][ti];
				Q.X[ti2][Ba[i]] = Q.X[ti][Ba[i]];

				Q.X[Ba[i]][bi2] = Q.X[Ba[i]][bi];
				Q.X[bi2][Ba[i]] = Q.X[bi][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][ti2] = Q.X[NBa[i]][ti];
				Q.X[ti2][NBa[i]] = Q.X[ti][NBa[i]];

				Q.X[NBa[i]][bi2] = Q.X[NBa[i]][bi];
				Q.X[bi2][NBa[i]] = Q.X[bi][NBa[i]];
			}
			Q.X[ti2][ti].setZero();
			Q.X[ti][ti2].setZero();
			Q.X[bi2][bi].setZero();
			Q.X[bi][bi2].setZero();
			Q.X[ti2][ti2] = Q.X[ti][ti];
			Q.X[bi2][bi2] = Q.X[bi][bi];
			Q.X[ti2][bi2] = Hi->q;
			Q.X[bi2][ti2] = Hi->q;

			ll bi2_idx = Ba.size() - 1;
			ll ti2_idx = NBa.size() - 1;
			ll bi_idx;
			ll ti_idx;
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == bi)
				{
					bi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == ti)
				{
					ti_idx = i;
					break;
				}
			}

			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in Hi\B*
				if (Hi_elements.find(NBa[i]) == Hi_elements.end())
				{
					C.X[bi2_idx][i] = C.X[ti_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (Hi_elements.find(Ba[i]) != Hi_elements.end())
				{
					C.X[i][ti2_idx] = C.X[i][bi_idx];
				}
			}
			p[bi2] = p[ti] + Q.X[bi2][ti];
			p[ti2] = p[bi] - Q.X[bi2][ti];

			M2.push_back(pr(ti_idx, bi_idx));
			M2.push_back(pr(bi2_idx, ti2_idx));
		}

		Pivoting_around_M(C, M2, Ba, NBa, Ba_sub, NBa_sub, inBase);
	}

	for (ll I = 0; I < Blossoms_in_MaximalNormalBlossom.size(); I++)
	{
		//node *Hi = NormalBlossoms[I];
		node *Hi = Blossoms_in_MaximalNormalBlossom[I];

		/*
		for (ll i = 0; i < Hi->ordering.size(); i++)
		{
			cout << Hi->ordering[i] << " ";
		}
		cout << endl;
		*/

		ll former_tip = Hi->tip;
		ll former_bud = Hi->bud;

		remove_vertex(Hi->bud, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
		remove_vertex(Hi->tip, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);

		Hi->tip = Hi->ordering[Hi->ordering.size() - 2];
		Hi->bud = Hi->ordering[Hi->ordering.size() - 1];
		Hi->ordering.erase(Hi->ordering.end() - 1);

		ll idx;
		for (ll j = 0; j < Hi->ordering.size(); j++)
		{
			if (Hi->ordering[j] == former_tip)
			{
				idx = j;
				break;
			}
		}

		Hi->ordering.erase(Hi->ordering.begin() + idx);

		vector<node *> Hi_ancestors = ancestors(Bl, Hi);
		for (ll j = 0; j < Hi_ancestors.size(); j++)
		{
			//cout << "4-" << j << endl;
			node *Hj = Hi_ancestors[j];
			/*
			for (ll k = 0; k < Hj->ordering.size(); k++)
			{
				cout << Hj->ordering[k] << " ";
			}
			cout << endl;
			*/

			ll idx;
			for (ll k = 0; k < Hj->ordering.size(); k++)
			{

				if (Hj->ordering[k] == former_tip)
				{
					idx = k;
					break;
				}
			}

			//cout << idx << endl;
			//cout << Hj->ordering.size() << endl;

			auto itr = Hj->ordering.erase(Hj->ordering.begin() + idx);

			for (ll k = 0; k < Hj->ordering.size(); k++)
			{
				if (Hj->ordering[k] == former_bud)
				{
					idx = k;
					break;
				}
			}

			itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
		}
	}
	/*cout<<"all blossom"<<endl;
	for(ll i=0;i<L.size();i++)
	{
		node *Hi=L[i];
		cout<<"key:"<<Hi->key<<endl;
		cout<<"tip:"<<Hi->tip<<" bud:"<<Hi->bud<<endl;
		cout<<"ordering:";
		for(ll j=0;j<Hi->ordering.size();j++)
		{
			cout<<Hi->ordering[j]<<" ";
		}
		cout<<endl;
	} */

	//-----step5-----
	//cout << "-----step5-----" << endl;
	vector<pr> new_dummy_line;
	for (ll I = 0; I < Blossoms_in_MaximalNormalBlossom.size(); I++)
	{
		//node *Hi = NormalBlossoms[I];
		node *Hi = Blossoms_in_MaximalNormalBlossom[I];

		ll ti2 = new_vertex(K, P, C, Q, order, rho, inBase, p);
		ll bi2 = new_vertex(K, P, C, Q, order, rho, inBase, p);

		vector<node *> Hi_ancestors = ancestors(Bl, Hi);

		Hi->ordering.push_back(ti2);
		new_dummy_line.push_back(pr(ti2, bi2));
		for (ll j = 0; j < Hi_ancestors.size(); j++)
		{
			node *Hj = Hi_ancestors[j];
			Hj->ordering.push_back(bi2);
			Hj->ordering.push_back(ti2);
		}
		set<ll> Hi_elements;
		for (ll j = 0; j < Hi->ordering.size(); j++)
		{
			Hi_elements.insert(Hi->ordering[j]);
		}
		ll ti = Hi->tip;
		ll bi = Hi->bud;

		vector<pr> M2;
		if (inBase[bi] && !inBase[ti])
		{
			Ba.push_back(ti2);
			NBa.push_back(bi2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[ti2] = true;
			inBase[bi2] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				Q.X[Ba[i]][ti2] = Q.X[Ba[i]][ti];
				Q.X[ti2][Ba[i]] = Q.X[ti][Ba[i]];

				Q.X[Ba[i]][bi2] = Q.X[Ba[i]][bi];
				Q.X[bi2][Ba[i]] = Q.X[bi][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][ti2] = Q.X[NBa[i]][ti];
				Q.X[ti2][NBa[i]] = Q.X[ti][NBa[i]];

				Q.X[NBa[i]][bi2] = Q.X[NBa[i]][bi];
				Q.X[bi2][NBa[i]] = Q.X[bi][NBa[i]];
			}
			Q.X[ti2][ti].setZero();
			Q.X[ti][ti2].setZero();
			Q.X[bi2][bi].setZero();
			Q.X[bi][bi2].setZero();
			Q.X[ti2][ti2] = Q.X[ti][ti];
			Q.X[bi2][bi2] = Q.X[bi][bi];
			Q.X[ti2][bi2] = Hi->q;
			Q.X[bi2][ti2] = Hi->q;

			ll ti2_idx = Ba.size() - 1;
			ll bi2_idx = NBa.size() - 1;
			ll bi_idx;
			ll ti_idx;
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == bi)
				{
					bi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == ti)
				{
					ti_idx = i;
					break;
				}
			}
			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in Hi\B*
				if (Hi_elements.find(NBa[i]) != Hi_elements.end())
				{
					C.X[ti2_idx][i] = C.X[bi_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (Hi_elements.find(Ba[i]) == Hi_elements.end())
				{
					C.X[i][bi2_idx] = C.X[i][ti_idx];
				}
			}
			p[bi2] = p[ti] - Q.X[bi2][ti];
			p[ti2] = p[bi] + Q.X[bi2][ti];

			M2.push_back(pr(bi_idx, ti_idx));
			M2.push_back(pr(ti2_idx, bi2_idx));
		}
		else if (!inBase[bi] && inBase[ti])
		{
			Ba.push_back(bi2);
			NBa.push_back(ti2);

			C.X.resize(C.row + 1);
			C.row++;
			for (ll i = 0; i < C.row; i++)
			{
				C.X[i].resize(C.col + 1);
			}
			C.col++;

			inBase[bi2] = true;
			inBase[ti2] = false;

			for (ll i = 0; i < Ba.size(); i++)
			{
				Q.X[Ba[i]][ti2] = Q.X[Ba[i]][ti];
				Q.X[ti2][Ba[i]] = Q.X[ti][Ba[i]];

				Q.X[Ba[i]][bi2] = Q.X[Ba[i]][bi];
				Q.X[bi2][Ba[i]] = Q.X[bi][Ba[i]];
			}
			for (ll i = 0; i < NBa.size(); i++)
			{
				Q.X[NBa[i]][ti2] = Q.X[NBa[i]][ti];
				Q.X[ti2][NBa[i]] = Q.X[ti][NBa[i]];

				Q.X[NBa[i]][bi2] = Q.X[NBa[i]][bi];
				Q.X[bi2][NBa[i]] = Q.X[bi][NBa[i]];
			}
			Q.X[ti2][ti].setZero();
			Q.X[ti][ti2].setZero();
			Q.X[bi2][bi].setZero();
			Q.X[bi][bi2].setZero();
			Q.X[ti2][ti2] = Q.X[ti][ti];
			Q.X[bi2][bi2] = Q.X[bi][bi];
			Q.X[ti2][bi2] = Hi->q;
			Q.X[bi2][ti2] = Hi->q;

			ll bi2_idx = Ba.size() - 1;
			ll ti2_idx = NBa.size() - 1;
			ll bi_idx;
			ll ti_idx;
			for (ll i = 0; i < NBa.size(); i++)
			{
				if (NBa[i] == bi)
				{
					bi_idx = i;
					break;
				}
			}
			for (ll i = 0; i < Ba.size(); i++)
			{
				if (Ba[i] == ti)
				{
					ti_idx = i;
					break;
				}
			}

			for (ll i = 0; i < C.col; i++)
			{
				//if NBa[i] in Hi\B*
				if (Hi_elements.find(NBa[i]) == Hi_elements.end())
				{
					C.X[bi2_idx][i] = C.X[ti_idx][i];
				}
			}
			for (ll i = 0; i < C.row; i++)
			{
				if (Hi_elements.find(Ba[i]) != Hi_elements.end())
				{
					C.X[i][ti2_idx] = C.X[i][bi_idx];
				}
			}
			p[bi2] = p[ti] + Q.X[bi2][ti];
			p[ti2] = p[bi] - Q.X[bi2][ti];

			M2.push_back(pr(ti_idx, bi_idx));
			M2.push_back(pr(bi2_idx, ti2_idx));
		}
		Pivoting_around_M(C, M2, Ba, NBa, Ba_sub, NBa_sub, inBase);
	}

	//cout << "-----step5 remove-----" << endl;
	for (ll I = 0; I < Blossoms_in_MaximalNormalBlossom.size(); I++)
	{
		//node *Hi = NormalBlossoms[I];
		node *Hi = Blossoms_in_MaximalNormalBlossom[I];

		ll former_tip = Hi->tip;
		ll former_bud = Hi->bud;
		remove_vertex(Hi->bud, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
		remove_vertex(Hi->tip, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);

		Hi->tip = new_dummy_line[I].first;
		Hi->bud = new_dummy_line[I].second;

		ll idx;
		for (ll j = 0; j < Hi->ordering.size(); j++)
		{
			if (Hi->ordering[j] == former_tip)
			{
				idx = j;
				break;
			}
		}

		Hi->ordering.erase(Hi->ordering.begin() + idx);

		vector<node *> Hi_ancestors = ancestors(Bl, Hi);
		for (ll j = 0; j < Hi_ancestors.size(); j++)
		{
			//cout << "3-" << j << endl;
			node *Hj = Hi_ancestors[j];
			ll idx;
			for (ll k = 0; k < Hj->ordering.size(); k++)
			{
				if (Hj->ordering[k] == former_tip)
				{
					idx = k;
					break;
				}
			}

			auto itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
			for (ll k = 0; k < Hj->ordering.size(); k++)
			{
				if (Hj->ordering[k] == former_bud)
				{
					idx = k;
					break;
				}
			}
			itr = Hj->ordering.erase(Hj->ordering.begin() + idx);
		}
	}

	/* cout<<"all blossom"<<endl;
	for(ll i=0;i<L.size();i++)
	{
		node *Hi=L[i];
		cout<<"key:"<<Hi->key<<endl;
		cout<<"tip:"<<Hi->tip<<" bud:"<<Hi->bud<<endl;
		cout<<"ordering:";
		for(ll j=0;j<Hi->ordering.size();j++)
		{
			cout<<Hi->ordering[j]<<" ";
		}
		cout<<endl;
	} */
	//cout << "C:";
	//C.output_matrix();
	//cout << endl;
	/*
	for (ll i = 0; i < Ba.size(); i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			ll v = Ba[i];
			ll u = NBa[j];

			ll v_idx, u_idx;
			for (ll k = 0; k < Ba.size(); k++)
			{
				if (Ba[k] == v)
				{
					v_idx = k;
					break;
				}
			}
			for (ll k = 0; k < NBa.size(); k++)
			{
				if (NBa[k] == u)
				{
					u_idx = k;
					break;
				}
			}
			if (!C.X[v_idx][u_idx].isZero())
			{
				if (p[u] - p[v] == Q.X[u][v])
				{
					cout << "tight:" << v << " " << u << endl;
				}
				else
				{
					cout << "no tight:" << v << " " << u << endl;
				}
			}
		}
	}
	*/
	//cout << "-----------Augment end-------" << endl;
}

void Search_in_Blossom(Matrix &C, node *H, Tree &Bl, vector<Field> &p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &order, vector<ll> &rho, vector<ll> &K, vector<bool> &inBase, Matrix &Q, vector<vector<ll>> &P, ll N)
{
	/////
	/*
	cout << endl;
	cout << "-------Search_in_Blossom start--------" << endl;
	cout << "blossom:" << H->key << endl;
	cout << "elements:";
	for (ll i = 0; i < H->ordering.size(); i++)
	{
		cout << H->ordering[i] << " ";
	}
	cout << endl;
	*/

	/*
	cout << "p:" << endl;
	for (ll i = 0; i < p.size(); i++)
	{
		cout << i << " ";
		p[i].output();
		cout << endl;
	}
	*/

	/*
	for (ll i = 0; i < Ba.size(); i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			ll v = Ba[i];
			ll u = NBa[j];

			ll v_idx, u_idx;
			for (ll k = 0; k < Ba.size(); k++)
			{
				if (Ba[k] == v)
				{
					v_idx = k;
					break;
				}
			}
			for (ll k = 0; k < NBa.size(); k++)
			{
				if (NBa[k] == u)
				{
					u_idx = k;
					break;
				}
			}
			if (!C.X[v_idx][u_idx].isZero())
			{
				if (p[u] - p[v] == Q.X[u][v])
				{
					cout << "tight:" << v << " " << u << endl;
				}
				else
				{
					cout << "no tight:" << v << " " << u << endl;
					p[u].output();
					cout << " ";
					p[v].output();
					cout << " ";
					Q.X[u][v].output();
					cout << endl;
				}
			}
		}
	}
	*/

	//vector<node *> all_blossom = all_blossoms(Bl);
	/* for (ll i = 0; i < all_blossom.size(); i++)
	{
		node *H_ = all_blossom[i];
		cout << "key:" << H_->key << " "
			 << "elements:";
		for (ll j = 0; j < H_->ordering.size(); j++)
		{
			cout << H_->ordering[j] << " ";
		}
		cout << endl;
		cout << "q:";
		H_->q.output();
		cout << endl;
	} */

	//cout << "C:" << endl;
	//C.output_matrix();
	//cout << endl;

	set<ll> H_elements;

	for (ll i = 0; i < H->ordering.size(); i++)
	{
		H_elements.insert(H->ordering[i]);
	}
	H_elements.insert(H->bud);

	Matrix C_H;
	vector<ll> H_Ba;
	vector<ll> H_Ba_idx;
	vector<ll> H_NBa;
	vector<ll> H_NBa_idx;

	for (ll i = 0; i < Ba.size(); i++)
	{
		if (H_elements.find(Ba[i]) != H_elements.end())
		{
			H_Ba.push_back(Ba[i]);
			H_Ba_idx.push_back(i);
		}
	}
	for (ll i = 0; i < NBa.size(); i++)
	{
		if (H_elements.find(NBa[i]) != H_elements.end())
		{
			H_NBa.push_back(NBa[i]);
			H_NBa_idx.push_back(i);
		}
	}

	Tree H_tree;
	H_tree.num = Bl.num;
	H_tree.root = H;

	while (!H->q.isZero())
	{
		//cout<<"q(Hi):"; H->q.output(); cout<<endl;
		Search(C, inBase, Ba, NBa, H_Ba, H_NBa, Q, K, H_tree, Bl, p, order, rho, P, N);

		//update dual variables
		vector<node *> MaximalBlossoms = maximal_blossom(H_tree);
		vector<ll> RZ_plus;
		RZ_plus.push_back(H->tip);
		for (ll i = 0; i < MaximalBlossoms.size(); i++)
		{
			node *Hi = MaximalBlossoms[i];
			if (Hi->label == 1)
			{
				for (ll j = 0; j < Hi->ordering.size(); j++)
				{
					RZ_plus.push_back(Hi->ordering[j]);
				}
			}
		}

		vector<ll> V_inverse;
		V_inverse.resize(order.size());

		for (ll i = 0; i < Ba.size(); i++)
		{
			V_inverse[Ba[i]] = i;
		}
		for (ll i = 0; i < NBa.size(); i++)
		{
			V_inverse[NBa[i]] = i;
		}
		Field eps;
		int ope = 0;
		for (ll i = 0; i < RZ_plus.size(); i++)
		{
			ll u = RZ_plus[i];
			for (ll j = 0; j < RZ_plus.size(); j++)
			{
				ll v = RZ_plus[j];
				if (inBase[u] && !inBase[v])
				{
					if (!C.X[V_inverse[u]][V_inverse[v]].isZero())
					{
						if (!judge_K(u, v, K))
						{
							if (ope == 0)
							{
								eps = p[v] - p[u] - Q.X[u][v];
								ope = 1;
							}
							else
							{
								eps = min(eps, (p[v] - p[u] - Q.X[u][v]));
							}
						}
					}
				}
				else if (!inBase[u] && inBase[v])
				{
					if (!C.X[V_inverse[v]][V_inverse[u]].isZero() && !judge_K(u, v, K))
					{
						if (ope == 0)
						{
							eps = p[u] - p[v] - Q.X[u][v];
							ope = 1;
						}
						else
						{
							eps = min(eps, p[u] - p[v] - Q.X[u][v]);
						}
					}
				}
			}
		}
		Field eps_;
		if (ope == 1)
		{
			eps_ = min(eps.divideByTwo(), H->q);
		}
		else
		{
			eps_ = H->q;
		}

		////
		/*
		cout << "eps_:";
		eps_.output();
		cout << endl; 
		*/
		/////

		if (inBase[H->tip])
		{
			p[H->tip] += eps_;
		}
		else
		{
			p[H->tip] -= eps_;
		}
		H->q -= eps_;

		for (ll i = 0; i < H->ordering.size(); i++)
		{
			ll v = H->ordering[i];
			for (ll j = 0; j < Ba.size(); j++)
			{
				Q.X[v][Ba[j]] -= eps_;
				Q.X[Ba[j]][v] -= eps_;
			}
			for (ll j = 0; j < NBa.size(); j++)
			{
				Q.X[v][NBa[j]] -= eps_;
				Q.X[NBa[j]][v] -= eps_;
			}
		}
		for (ll i = 0; i < H->ordering.size(); i++)
		{
			for (ll j = 0; j < H->ordering.size(); j++)
			{
				ll u = H->ordering[i];
				ll v = H->ordering[j];
				if (u != v)
				{
					Q.X[u][v] += eps_;
					Q.X[v][u] += eps_;
				}
				else
				{
					Q.X[u][u] += eps_;
				}
			}
		}

		for (ll i = 0; i < MaximalBlossoms.size(); i++)
		{
			node *Hi = MaximalBlossoms[i];
			Hi->q += eps_;

			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				ll v = Hi->ordering[j];
				for (ll k = 0; k < Ba.size(); k++)
				{
					Q.X[v][Ba[k]] += eps_;
					Q.X[Ba[k]][v] += eps_;
				}
				for (ll k = 0; k < NBa.size(); k++)
				{
					Q.X[v][NBa[k]] += eps_;
					Q.X[NBa[k]][v] += eps_;
				}
			}
			for (ll j = 0; j < Hi->ordering.size(); j++)
			{
				for (ll k = 0; k < Hi->ordering.size(); k++)
				{
					ll u = Hi->ordering[j];
					ll v = Hi->ordering[k];
					if (u != v)
					{
						Q.X[u][v] = Q.X[u][v] - eps_;
						Q.X[v][u] = Q.X[v][u] - eps_;
					}
					else
					{
						Q.X[u][u] -= eps_;
					}
				}
			}

			if (inBase[Hi->bud])
			{
				p[Hi->bud] -= eps_;
			}
			else
			{
				p[Hi->bud] += eps_;
			}
		}
		/*
		vector<node*> blossoms_in_H=all_blossoms(H_tree);
		for(ll i=0;i<blossoms_in_H.size();i++){
			node *Hi=blossoms_in_H[i];
			if(Hi->q.isZero())
			{
				Expand(C, Hi, Ba, NBa, H_Ba, H_NBa, inBase, order, Bl);
			}
		}
		*/
	}

	Expand(C, H, Ba, NBa, H_Ba, H_NBa, inBase, order, Bl);

	/* vector<node *> allblo = all_blossoms(Bl);
	cout << "BLOSSOM:";
	for (ll i = 0; i < allblo.size(); i++)
	{
		node *Hi = allblo[i];
		cout << Hi->key << " ";
	}
	cout << endl;
	/*/
	//cout << "-------Search_in_Blossom end----------" << endl;
	//cout << endl;
}

vector<ll> matroid_parity(Matrix &A, Field w[])
{
	//Aから行列Cをつくる
	Matrix A_b(A.row, A.row);
	Matrix A_nb(A.row, A.col - A.row);

	vector<Field> p;
	p.resize(A.col);
	for (ll i = 0; i < A.col; i++)
	{
		p[i] = w[i / 2].divideByTwo();
	}

	//cout << "Base_Greedy start" << endl;
	vector<ll> Ba = Base_Greedy(A, p);

	cout << "Base:";
	for (ll i = 0; i < Ba.size(); i++)
	{
		cout << Ba[i] << " ";
	}
	cout << endl;
	//cout << "Base Greedy end" << endl;
	if (overflow)
	{
		return {};
	}
	if (Ba.size() != A.row)
	{
		//cout << "no base" << endl;
		return Ba;
	}
	//cout<<"0"<<endl;

	vector<bool> inBase;
	inBase.resize(A.col, false);
	for (ll i = 0; i < Ba.size(); i++)
	{
		inBase[Ba[i]] = true;
	}
	vector<ll> NBa;
	for (ll i = 0; i < A.col; i++)
	{
		if (inBase[i] == false)
		{
			NBa.push_back(i);
		}
	}
	/*
	cout << "NBa:";
	for (ll i = 0; i < NBa.size(); i++)
	{
		cout << NBa[i] << " ";
	}
	cout << endl;
	*/

	for (ll i = 0; i < A.row; i++)
	{
		for (ll j = 0; j < Ba.size(); j++)
		{
			A_b.X[i][j] = A.X[i][Ba[j]];
		}
	}
	for (ll i = 0; i < A.row; i++)
	{
		for (ll j = 0; j < NBa.size(); j++)
		{
			A_nb.X[i][j] = A.X[i][NBa[j]];
		}
	}

	//cout<<"1"<<endl;
	Matrix C = matrix_inverse(A_b) * A_nb;
	//cout<<"11"<<endl;
	//cout<<"C:"<<endl;
	//C.output_matrix();
	/* 	if (overflow)
	{
		overflow = true;
		return {};
	} */
	num_type C_mx = C.maximum_number();
	//cout << "the maximum absolute number in elements of C:";
	//cout << C_mx << endl;

	maximum_absolute_num = max(maximum_absolute_num, C_mx);
	//cout<<"111"<<endl;

	//matrix_inverse(A_b).output_matrix();

	vertex_number = A.col;

	vector<ll> K;
	K.resize(A.col);
	for (ll i = 0; i < A.col; i++)
	{
		K[i] = -i;
	}
	Matrix Q(A.col, A.col); //Q[i][i] : sum of the value q over H where i is included
	for (ll i = 0; i < Q.row; i++)
	{
		for (ll j = 0; j < Q.col; j++)
		{
			Q.X[i][j].setZero();
		}
	}

	Tree Blossom(A.col);

	vector<ll> order;
	vector<ll> rho;
	vector<vector<ll>> P;

	/* bool flag = SourceLine_vector(A.col, Ba, inBase);
	if (flag == true)
	{
		cout << "source line" << endl;
	}
	else
	{
		cout << "no source line" << endl;
	}

	cout << "step2 start" << endl; */
	while (SourceLine_vector(A.col, Ba, inBase))
	{
		//cout<<"2"<<endl;
		//cout << "Search start" << endl;

		vector<ll> Ba_sub = Ba;
		vector<ll> NBa_sub = NBa;

		vector<ll> AugmentingPath = Search(C, inBase, Ba, NBa, Ba_sub, NBa_sub, Q, K, Blossom, Blossom, p, order, rho, P, A.col);
		//cout << "Search end" << endl;
		if (overflow)
		{
			return {};
		}

		if (AugmentingPath.empty())
		{
			//cout << "DualUpdate start" << endl;
			if (DualUpdate(C, Blossom, p, Ba, NBa, order, rho, K, inBase, Q))
			{
				/*
				vector<node *> AllBlossoms = all_blossoms(Blossom);
				for (ll i = 0; i < AllBlossoms.size(); i++)
				{
					node *Hi = AllBlossoms[i];
					if (Hi->q.isZero())
					{
						//cout << "Expand H" << Hi->key << endl;
						Expand(C, Hi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Blossom);
						if (overflow)
						{
							return {};
						}
					}
				}
				*/

				/////変えた
				/*
				while(1){
					vector<node*>MaximalBlossom=maximal_blossom(Blossom);
					if(MaximalBlossom.size()==0){
						break;
					}
					for(ll i=0;i<MaximalBlossom.size();i++)
					{
						node *Hi=MaximalBlossom[i];
						if(Hi->q.isZero())
						{
							Expand(C, Hi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Blossom);
						}
					}
				}
				*/
				vector<node *> AllBlossoms = all_blossoms(Blossom);
				//大きい花から小さい花へ
				for (ll i = 0; i < AllBlossoms.size(); i++)
				{
					node *Hi = AllBlossoms[i];
					if (Hi->q.isZero() && isMaximal(Blossom, Hi))
					{
						//cout << "Expand H" << Hi->key << endl;
						Expand(C, Hi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Blossom);
						if (overflow)
						{
							return {};
						}
					}
				}
				/////
			}
			else
			{
				return {};
			}
		}
		else
		{
			/* cout << "DualUpdate end" << endl;
			cout << "Augment start" << endl; */
			Augment(AugmentingPath, Blossom, C, Ba, NBa, Ba_sub, NBa_sub, P, K, order, rho, inBase, Q, p);
			if (overflow)
			{
				return {};
			}
			//cout << "Augment end" << endl;
			vector<node *> AllBlossoms = all_blossoms(Blossom);

			//cout << "Search in each blossom start" << endl;
			vector<node *> Blossoms_in_MaximalNormalBlossom = blossoms_in_MaximalNormalBlossom(Blossom);
			/*for (ll i = AllBlossoms.size() - 1; i >= 0; i--)
			{
				node *Hi = AllBlossoms[i];
				if (Hi->normal)
				{
					Search_in_Blossom(C, Hi, Blossom, p, Ba, NBa, order, rho, K, inBase, Q, P, A.col);
				}
			}
			*/
			for (ll i = Blossoms_in_MaximalNormalBlossom.size() - 1; i >= 0; i--)
			{
				node *Hi = Blossoms_in_MaximalNormalBlossom[i];
				Search_in_Blossom(C, Hi, Blossom, p, Ba, NBa, order, rho, K, inBase, Q, P, A.col);
				if (overflow)
				{
					return {};
				}
			}
			//cout << "Search_in_Blossom end" << endl;
		}
	}

	//cout << "main routine end" << endl;
	return Ba;
}

pair<Matrix, vector<RationalNumber>> create_matrix_for_spath(vector<pair<pair<int, int>, RationalNumber>> edges, int path_num, int vertex_num, int edge_num, vector<vector<int>> A, int terminal_num)
{
	vector<RationalNumber> w;
	int additional_terminal_num = terminal_num - 2 * path_num;

	Matrix incidence_matrix(vertex_num + additional_terminal_num + 2, edge_num + terminal_num * additional_terminal_num + (vertex_num - terminal_num) + 1);

	for (int i = 0; i < edges.size(); i++)
	{
		incidence_matrix.X[edges[i].first.first][i] = RationalNumber(1, 1);
		incidence_matrix.X[edges[i].first.second][i] = RationalNumber(-1, 1);
		w.push_back(edges[i].second);
	}

	int cnt = 0;
	A.resize(A.size() + 1);

	for (int i = 0; i < additional_terminal_num; i++)
	{
		for (int j = 0; j < terminal_num; j++)
		{
			for (int l = 0; l < A[j].size(); l++)
			{
				//terminals.insert(A[j][l]);
				incidence_matrix.X[A[j][l]][edge_num + cnt] = RationalNumber(1, 1);
				incidence_matrix.X[vertex_num + i][edge_num + cnt] = -RationalNumber(-1, 1);
				cnt++;
				w.push_back(RationalNumber(0, 1));
			}
		}
		A[A.size() - 1].push_back(vertex_num + i);
		//terminals.insert(n+i);
	}

	incidence_matrix.X[vertex_num + additional_terminal_num][edge_num + cnt + 1] = RationalNumber(1, 1);
	incidence_matrix.X[vertex_num + additional_terminal_num + 1][edge_num + cnt + 1] = RationalNumber(-1, 1);

	A.resize(A.size() + 1);
	A[A.size() - 1].push_back(vertex_num + additional_terminal_num);
	w.push_back(RationalNumber(0, 1));
	A.resize(A.size() + 1);
	A[A.size() - 1].push_back(vertex_num + additional_terminal_num + 1);
	w.push_back(RationalNumber(0, 1));
	//terminals.insert(n+ad);
	//terminals.insert(n+ad+1);

	Matrix E2(2, 2);
	E2.X[0][0].setOne();
	E2.X[0][1].setZero();
	E2.X[1][0].setZero();
	E2.X[1][1].setOne();
	incidence_matrix = kronecker(incidence_matrix, E2);

	//消していく
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			int terminal = A[i][j];
			for (int l = 0; l < incidence_matrix.col; l++)
			{
				if (incidence_matrix.X[terminal * 2][l] == RationalNumber(1, 1))
				{
					incidence_matrix.X[terminal * 2][l].setZero();
					incidence_matrix.X[terminal * 2 + 1][l] = RationalNumber(-(i + 1), 1);
				}
				else if (incidence_matrix.X[terminal * 2][l] == RationalNumber(-1, 1))
				{
					incidence_matrix.X[terminal * 2][l].setZero();
					incidence_matrix.X[terminal * 2 + 1][l] = RationalNumber(i + 1, 1);
				}
			}
		}
	}

	return make_pair(incidence_matrix, w);

	/*
    for(int i=0;i<incidence_matrix.col;i++){
    	int node_num=incidence_matr
        if(i<edges.size()){
            if(terminals.find(edges[i].first)==terminals.end()){
                if(terminals.find(edges[i].second)!=terminals.end()){

                }
            }
            else{
                if(terminals.find(edges[i].second)==terminals.end()){

                }else{

                }
            }
        }
        else{

        }
    }
    */
}

/*
int main()
{
	vector<pair<pair<int,int>, RationalNumber> > edges;
	RationalNumber r=RationalNumber(1,1);
	edges.push_back(make_pair(make_pair(0,2),r));
	edges.push_back(make_pair(make_pair(1,2),r));
	edges.push_back(make_pair(make_pair(0,3),r));
	edges.push_back(make_pair(make_pair(1,3),r));
	edges.push_back(make_pair(make_pair(1,4),r));
	edges.push_back(make_pair(make_pair(2,3),r));
	edges.push_back(make_pair(make_pair(3,4),r));
	edges.push_back(make_pair(make_pair(0,4),r));
	int k=2;
	int vertex_num=5;
	int edge_num=8;
	int terminal_num=3;
	vector<vector<int> > A;
	A.resize(3);
	A[0].push_back(0);
	A[1].push_back(1);
	A[2].push_back(2);

	pair<Matrix,vector<RationalNumber> > p=create_matrix_for_spath(edges,k,vertex_num,edge_num,A,terminal_num);
	RationalNumber* w=new RationalNumber[p.first.col/2];
	cout<<"matrix"<<endl;
	p.first.output_matrix();
	cout<<"w:";
	for(int i=0;i<p.second.size();i++){
		w[i]=p.second[i];
		w[i].output();cout<<" ";
		}
		cout<<endl;
	vector<ll> minimum_weight_parity_base=matroid_parity(p.first,w);
	cout<<minimum_weight_parity_base.size()<<endl;
	for(int i=0;i<minimum_weight_parity_base.size();i++){
		cout<<minimum_weight_parity_base[i]<<" ";
	}

	return 0;
}
*/

int main()
{

	int mode;
	cout << "ランダムなら０，接続行列なら1，接続行列の変種なら２を入力" << endl;
	cin >> mode;
	//ll element = 1;
	//ll weight = 10;
	//double nonzero_prob = 0.5;
	//vector<int> row_size = {10,20,50,100,200};
	//vector<int> row_size={200};
	//vector<int> row_size = {200,500};

	vector<int> elements;
	vector<int> row_size;
	const ll weight = 10;
	vector<string> nonzero_probs;

	if (mode == 0)
	{
		elements = {1, 2, 5, 10};
		row_size = {10};
		nonzero_probs = {"0.1", "0.2", "0.5"};
	}
	if (mode == 1)
	{
		elements = {1};
		row_size = {10, 20, 50, 100, 200};
		nonzero_probs = {"N"};
	}
	if (mode == 2)
	{
		elements = {3, 5, 10, 20, 30};
		row_size = {10, 20, 50, 100};
		nonzero_probs = {"N"};
	}
	if (mode == 4)
	{
		elements = {1};
		row_size = {10, 20, 50, 100, 200, 500, 1000, 2000};
		nonzero_probs = {"N"};
	}

	for (ll I = 0; I < row_size.size(); I++)
	{
		for (ll nonzero_idx = 0; nonzero_idx != nonzero_probs.size(); nonzero_idx++)
		{

			string nonzero_prob = nonzero_probs[nonzero_idx];
			for (ll element_idx = 0; element_idx != elements.size(); element_idx++)
			{
				ll element = elements[element_idx];
				for (ll file_num = 1; file_num <= 10; file_num++)
				{

					Blossom_num = 0;
					Graft_num = 0;
					Augment_num = 0;
					Search_num = 0;
					DualUpdate_num = 0;

					maximum_absolute_num = 0;
					overflow = false;

					string filepath;
					ll rowsize = row_size[I];

					if (mode == 0)
					{
						filepath = "./sotsuron_data/data_";
						filepath += to_string(weight) + "_" + to_string(element) + "_" + nonzero_prob;
						filepath += "/data_" + to_string(rowsize) + "_" + to_string(rowsize * 2) + "_" + to_string(weight) + "_" + to_string(element) + "_" + nonzero_prob + "_" + to_string(file_num) + ".txt";
						cout << filepath << endl;
					}
					if (mode == 1)
					{
						filepath = "./sotsuron_data/incidence_un_3/incidence_matrix_";
						filepath += to_string(rowsize) + "_" + to_string(rowsize * 3) + "_" + to_string(file_num) + ".txt";
						cout << filepath << endl;
					}
					if (mode == 2)
					{
						filepath = "./sotsuron_data/spath_";
						filepath += to_string(element) + "_0.5/spath_";
						filepath += to_string(element) + "_" + to_string(rowsize) + "_" + to_string(rowsize * 3) + "_" + to_string(file_num) + "_" + "0.5" + ".txt";
					}
					if (mode == 4)
					{
						filepath = "./sotsuron_data/yuukou_inc/yuukou_incidence_matrix_";
						filepath += to_string(rowsize) + "_" + to_string(rowsize * 3) + "_" + to_string(file_num) + ".txt";
						cout << filepath << endl;
					}

					ifstream ifs;
					ifs.open(filepath);

					if (ifs.fail())
					{
						cout << "Failed to open file." << endl;
						return -1;
					}

					string x, y;
					ifs >> x >> y;
					//cout << "x" << x << y << endl;
					Matrix A(stoll(x), stoll(y));
					//cout << A.row << " " << A.col << endl;
					A.matrix_setZero();
					//A.output_matrix();
					ll n = A.col / 2;
					Field *w = new Field[n];
					for (ll i = 0; i < n; i++)
					{
						w[i].setZero();
					}

					for (ll i = 0; i < A.col; i++)
					{
						string str;
						ifs >> str;
						ll num = stoll(str);

						A.X[0][i].input_ll(num);
						//A.X[0][i].output();
						//cout<<endl;
					}
					ll cnt = 0;
					while (!ifs.eof())
					{
						for (ll i = 0; i < A.col; i++)
						{
							string str;
							ifs >> str;
							ll num = stoll(str);
							A.X[cnt + 1][i].input_ll(num);
							//A.X[cnt+1][i].output();cout<<" ";
						}
						//cout<<endl;
						cnt++;
						//cout<<cnt<<endl;
						if (cnt == A.row - 1)
						{
							//cout<<"b"<<endl;
							for (ll i = 0; i < A.col / 2; i++)
							{
								string str;
								ifs >> str;
								ll num = stoll(str);
								w[i].input_ll(num);
								//cout<<"b-"<<i<<endl;
							}
							break;
						}
					}

					//cout<<"aaa"<<endl;
					//A.output_matrix();
					/*
			cout<<"w:";
			for(ll i=0;i<n;i++){
				w[i].output();cout<<" ";
			}
			cout<<endl;
			*/

					//cout << "行列の入力：" << endl;
					//Matrix A;
					//A.input_matrix_int();
					//cout << "a" << endl;

					/* ll n = A.col / 2;
		Field *w = new Field[n];
		for (ll i = 0; i < n; i++)
		{
			w[i].setZero();
		} */

					/*
		vector<ll> minimum_weight_parity_base = matroid_parity(A, w);
		for (ll i = 0; i < minimum_weight_parity_base.size(); i++)
		{
			cout << minimum_weight_parity_base[i] << " ";
		}
		if (minimum_weight_parity_base.size() == 0)
		{
			cout << "there is no parity base" << endl;
		}
		cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	}while(next_permutation(vv.begin(),vv.end()));
	*/

					//cout << "ラインの重み入力" << endl;
					/* for (ll i = 0; i < A.col / 2; i++)
				{
			//Field a;
			//a.input();
			//w[i] = a;
			ll a;
			cin >> a;
			w[i] = RationalNumber(a, 1);
						} */
					//cout << "b" << endl;

					clock_t start = clock();
					vector<ll> minimum_weight_parity_base = matroid_parity(A, w);
					clock_t end = clock();

					if (mode == 0)
					{
						cout << "size of A:" << A.row << "*" << A.col << " w:" << weight << " gamma:" << element << " nonzero prob:" << nonzero_prob << endl;
					}
					if (mode == 1)
					{
						cout << "size of A:" << A.row << "*" << A.col << " w:" << weight << endl;
					}
					if (mode == 2)
					{
						cout << "size of A:" << A.row << "*" << A.col << " w:" << weight << endl;
						cout << "k:" << element << endl;
					}
					if (mode == 4)
					{
						cout << "size of A:" << A.row << "*" << A.col << " w:" << weight << endl;
					}

					if (minimum_weight_parity_base.size() != 0)
					{
						sort(minimum_weight_parity_base.begin(), minimum_weight_parity_base.end());
					}

					cout << "the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment" << endl;
					cout << Search_num << " " << Blossom_num << " " << Graft_num << " " << DualUpdate_num << " " << Augment_num << endl;

					cout << "maximum absolute number in C:" << maximum_absolute_num << endl;
					cout << "solution:";
					if (minimum_weight_parity_base.size() != 0 && minimum_weight_parity_base.size() < A.row && !overflow)
					{
						cout << "A is not row full rank matrix" << endl;
						cout << "time:" << (double)(end - start) / CLOCKS_PER_SEC << "sec" << endl;
						cout << endl;
						continue;
					}

					Field opt_value;
					opt_value.setZero();
					for (ll i = 0; i < minimum_weight_parity_base.size(); i++)
					{
						if (minimum_weight_parity_base[i] < A.col)
						{
							cout << minimum_weight_parity_base[i] << " ";
							opt_value += w[minimum_weight_parity_base[i] / 2];
						}
					}
					cout << endl;
					if (minimum_weight_parity_base.size() != 0)
					{
						cout << "optimal value:";
						(opt_value.divideByTwo()).output();
						cout << endl;
					}
					if (minimum_weight_parity_base.size() == 0 && !overflow)
					{
						cout << "there is no parity base" << endl;
					}
					if (minimum_weight_parity_base.size() == 0 && overflow)
					{
						cout << "overflow" << endl;
					}

					cout << "time:" << (double)(end - start) / CLOCKS_PER_SEC << "sec" << endl;
					cout << endl;
				}

				cout << "---------------------------------------------" << endl;
				cout << endl;

				/* Matrix B(A.row, A.row);
	for (ll i = 0; i < A.row; i++)
	{
		for (ll j = 0; j < minimum_weight_parity_base.size(); j++)
		{
			B.X[i][j] = A.X[i][minimum_weight_parity_base[j]];
		}
	}
	cout << endl;
	matrix_inverse(B).output_matrix(); */
			}
			cout << "***********************************************************" << endl;
			cout << "***********************************************************" << endl;
			cout << endl;
		}
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << endl;
	}

	return 0;
}
