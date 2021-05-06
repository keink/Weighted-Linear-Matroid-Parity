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

num_type maximum_absolute_num = 0; // maximum absolute number that appear in a target matrix during executing an algorithm
bool overflow = false;

// Pivoting around a pair p
void Pivoting_around_p(Matrix &C, pr p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase)
{
	ll x = p.first;
	ll y = p.second; //x in B and y not in B

	for (ll i = 0; i < C.row; i++)
	{
		for (ll j = 0; j < C.col; j++)
		{
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

	maximum_absolute_num = max(maximum_absolute_num, C_mx);
}

// Pivoting around a matching M
void Pivoting_around_M(Matrix &C, vector<pr> &M, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase)
{
	vector<bool> flag_Ba;
	vector<bool> flag_NBa;
	flag_Ba.resize(C.row, false);
	flag_NBa.resize(C.col, false);
	vector<ll> A;
	vector<ll> B;

	for (ll i = 0; i < M.size(); i++)
	{
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

	Matrix alpha_inverse = matrix_inverse(alpha);

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
	}

	num_type C_mx = C.maximum_number();


	maximum_absolute_num = max(maximum_absolute_num, C_mx);

}

// Judges if {i, j} is a line or not
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

// Returns the mate of a given variable i
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

// Judges if the vertex i is a single vertex or not
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

// Judges if K(u) and K(v) matches each other
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

// Judges if there exists a source line in a given vector
// Note: N is the number of columns of A
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

// Judges if there exists a source line in a given set
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

// Gets the path from v to a source vertex
vector<ll> Path(ll v, vector<vector<ll>> &P)
{
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

	return path;
}

bool comp(const pr &a, const pr &b)
{
	return a.first > b.first;
}

// Adds new vertex to V
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

// Removes vertex from V
void remove_vertex(ll u, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, Matrix &C)
{
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
}

// Expands a blossom H
void Expand(Matrix &C, node *H, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, Tree &Bl)
{
	if (!H->normal)
	{
		tree_delete(Bl, H->key);
	}
	else
	{
		ll t = H->tip;
		ll b = H->bud;

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
		}

		Pivoting_around_p(C, p, Ba, NBa, Ba_sub, NBa_sub, inBase);

		// ↓ The following code is relevant  only when "Search in Blossom" is being executed. 
		vector<node *> H_ancestors = ancestors(Bl, H);
		for (ll i = 0; i < H_ancestors.size(); i++)
		{
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
}

// Judges if there exists an edge between u and v
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

		// Multiplies by a so that A[j][j]=1
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

	return ret;
}

// Gets a base of a given matrix A by greedy method
vector<ll> Base_Greedy(Matrix &A, vector<Field> &p)
{
	Matrix B = A;
	vector<pair<Field, ll>> M;
	for (ll i = 0; i < A.col; i++)
	{
		M.push_back(make_pair(p[i], i));
	}
	sort(M.begin(), M.end());

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

void Blossom(ll v, ll u, Matrix &C, Matrix &Q, vector<vector<ll>> &P, Tree &Bl, Tree &Bl2, vector<ll> &K, vector<ll> &order, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<Field> &p, vector<ll> &rho, queue<ll> &que, ll N)
{
	if (!overflow)
	{
		Blossom_num++;
	}

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
			// In this case, cur is a single vertex.
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

	ll r = -1;
	ll r_idx = -1;
	set<ll> children;
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		if (isSingle(path_v[i], K))
		{
			children.insert(path_v[i]); 
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

		// Wants to find r
		node *x = H->fchild;
		while (x != NULL)
		{
			//TODO: update comment
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
		else
		{
			children.insert(K[c]);
		}

		Field q_H;
		q_H.setZero();
		blossom_insert(Bl, blossom_number, children, q_H);
		blossom_number++;

		H = tree_search(Bl, blossom_number - 1 + Bl.num).first;

		// Wants to find r
		node *x = H->fchild;
		while (x != NULL)
		{
			//TODO:update comment
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
	H->routing.resize(vertex_number + 2);

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

		// ↓ The following code is relevant  only when "Search in Blossom" is being executed. 
		vector<node *> H_ancestors = ancestors(Bl2, H);
		for (ll i = 0; i < H_ancestors.size(); i++)
		{
			node *Hi = H_ancestors[i];
			Hi->ordering.push_back(H->tip);
			Hi->ordering.push_back(H->bud);
		}

		ll t = H->tip;
		ll b = H->bud;
		H->normal = true;

		for (ll i = 0; i < Ba.size(); i++)
		{
			// It's OK since q(H) will be set as 0
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

		// ↓ The following code is relevant  only when "Search in Blossom" is being executed. 
		if (Bl.root->key != -1)
		{
			node *Hj = Bl.root;
			ll tj = Hj->tip;

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
		}

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

			p[b] = p[r] - Q.X[r][b];
			p[t] = p[b];
		}

		Pivoting_around_p(C, pr(Ba.size() - 1, NBa.size() - 1), Ba, NBa, Ba_sub, NBa_sub, inBase);

		//-----step3-----
		if (!P[mate(g)].empty())
		{
			P[mate(g)][0] = t;
		}
		else
		{
			vector<node *> des = descendant(Bl, H);
		}

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

		// Labels t with P(t)=P(r)bt
		P[t].push_back(r);
		P[t].push_back(b);

		// Extends the ordering
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
	vector<pr> unlabeled;
	vector<pr> labeled;

	for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
	{
		if (order[*itr] < INF)
		{
			labeled.push_back(pr(order[*itr], *itr));
		}
	}

	// Sorts by values of order (descending order)
	sort(labeled.begin(), labeled.end(), comp);

	// Determins <H for labeled vertices
	for (ll i = labeled.size() - 1; i >= 0; i--)
	{
		H->ordering.push_back(labeled[i].second);
	}

	node *x = H->fchild;
	queue<pr> que_g;
	queue<pr> que_h;

// TODO: update comments
	//Step4の3，4番目の条件でlabelする頂点について
	//1,2番目の条件によりtiはすでにlabelされ，さらにRHi(x)もさだまっているから
	//P(x)=RHi(x)でよい
	//この後で必要なのはqueに入れる順番決め

	while (x != NULL)
	{
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
			int ope_h = 0;
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

	//unlabeledの頂点をP(x)でラベル
	//TODO: update comment
	for (ll i = d_idx + 1; i < path_u.size(); i++)
	{
		ll k = path_u[i];
		if (order[k] == INF && isSingle(k, K))
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

			unlabeled.push_back(pr(order[rho[k]], k));
		}
	}
	for (ll i = c_idx + 1; i < path_v.size(); i++)
	{
		ll k = path_v[i];
		if (order[k] == INF && isSingle(k, K))
		{
			if (i + 2 >= path_v.size())
			{
				if (u != r)
				{
					P[k].push_back(u);
				}
				else
				{ //normal line 一本のみの花をつくるとき
				//TODO: update comment
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

	if (!SourceLine_set(N, H_elements, inBase))
	{
		for (auto itr = H_elements.begin(); itr != H_elements.end(); itr++)
		{
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
}

void Graft(ll v, node *Hi, Matrix &C, Matrix &Q, vector<vector<ll>> &P, Tree &Bl, Tree &Bl2, vector<ll> &K, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<bool> &inBase, vector<ll> &order, vector<Field> &p, vector<ll> &rho, queue<ll> &que)
{
	if (!overflow)
	{
		Graft_num++;
	}
	//cout << "Graft start" << endl;
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
		// It's OK since q(H) will be set as 0
		Q.X[Ba[i]][t] = Q.X[Ba[i]][ti];
		Q.X[t][Ba[i]] = Q.X[ti][Ba[i]];
		// b is not in blossom cause H_i is a maximal blossom
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

	// Intoroduces b,t as well as Blossom step2
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
			// if NBa[i] in H\B*
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

		// Updates p
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
			// if NBa[i] in H\B*
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

		// Updates p
		p[b] = p[r] - Q.X[r][b];
		p[t] = p[b];
	}

	Pivoting_around_p(C, pr(Ba.size() - 1, NBa.size() - 1), Ba, NBa, Ba_sub, NBa_sub, inBase);

	// Labels t with P(t)=P(v)bt
	P[t].push_back(v);
	P[t].push_back(b);

	// Extends the ordering of < so that t is just after v
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
	Hi->label = 1;
	Hi->ordering[0] = t;
	// Defines RH(x) for each x in H\{bi}
	Hi->routing.resize(vertex_number);
	Hi->routing[t] = {t};

	Hi->bud = b;
	Hi->tip = t;
	K[b] = Hi->key;
	K[t] = Hi->key;

	//-----step4-----
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

	// ↓ The following code is relevant  only when "Search in Blossom" is being executed. 
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

	vector<ll> path;
	ll Num = inBase.size();

	queue<ll> que;

	// K[i] -i: single, number(0-|blossoms|) : maximal blossom
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

	// Initializes ordering with INF
	order.resize(Num, INF);
	ordering_number = 0;

	// Initializes K
	blossom_initialize(Bl, K);

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

	// unlabeled maximal source blossom
	vector<node *> MaximalBlossom = maximal_blossom(Bl);
	for (ll i = 0; i < MaximalBlossom.size(); i++)
	{

		node *x = MaximalBlossom[i];
		
		if (!x->normal)
		{
			// Labels maximal blossom with plus
			x->label = 1;
			for (ll j = 0; j < x->ordering.size(); j++)
			{
				P[x->ordering[j]].reserve(x->routing[j].size());
				P[x->ordering[j]] = x->routing[x->ordering[j]];
				
				que.push(x->ordering[j]);
				order[x->ordering[j]] = ordering_number;
				ordering_number++;
			}
		}
	}

	for (auto itr = source_lines.begin(); itr != source_lines.end(); itr++)
	{
		pr sourceline = *itr;
		ll x = sourceline.first;
		ll y = sourceline.second;

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
			if (p[y] - p[x] == Q.X[x][y])
			{
				Blossom(x, y, C, Q, P, Bl, Bl2, K, order, Ba, NBa, Ba_sub, NBa_sub, inBase, p, rho, que, N);
			}
		}
	}

	//-----step2-----
	while (!que.empty())
	{
		ll v = que.front();
		que.pop();

		ll v_idx;

		while (1)
		{
			ll u;
			ll mn_order = INF;
			ll mn_order_idx;

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
									if (mn_order > order[NBa[i]])
									{
										u = NBa[i];
										mn_order = order[NBa[i]];
										mn_order_idx = i;
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
									if (mn_order > order[Ba[i]])
									{
										u = Ba[i];
										mn_order = order[Ba[i]];
										mn_order_idx = i;
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
								unlabeled.push_back(Ba[i]);
							}
						}
					}
					j++;
				}
			}
		}

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
	//cout << "-----------Search end------------------" << endl;
	return {};
}

bool DualUpdate(Matrix &C, Tree &Bl, vector<Field> &p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &order, vector<ll> &rho, vector<ll> &K, vector<bool> &inBase, Matrix &Q)
{
	if (!overflow)
	{
		DualUpdate_num++;
	}
	//cout << "DualUpdate start" << endl;

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

	// Defines Y
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

	// Defines eps
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

	// Gets eps1
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

	// Gets eps2
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

	// Gets eps3
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

	// Gets eps4
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

	Field q;
	q.setZero();
	if (eps < q)
	{
		exit(EXIT_FAILURE);
	}

	// Updates dual variables
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
				}
				for (ll k = 0; k < NBa.size(); k++)
				{
					Q.X[v][NBa[k]] += eps;
					Q.X[NBa[k]][v] += eps;
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
						Q.X[u][v] = Q.X[u][v] - eps;
						Q.X[v][u] = Q.X[v][u] - eps;
					}
					else
					{
						// Substracts the extra value added for Q_uu
						Q.X[u][u] -= eps;
					}
				}
			}
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
		}
	}

	return true;
}

void Augment(vector<ll> &path, Tree &Bl, Matrix &C, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &Ba_sub, vector<ll> &NBa_sub, vector<vector<ll>> &P, vector<ll> &K, vector<ll> &order, vector<ll> &rho, vector<bool> &inBase, Matrix &Q, vector<Field> &p)
{
	if (!overflow)
	{
		Augment_num++;
	}
	//cout << "-------Augment start-----------" << endl;

	//-----step0-----
	vector<node *> L = all_blossoms(Bl);

	set<ll> Lp;
	vector<node *> NLp;
	vector<pair<node *, pair<ll, ll>>> Lp_positive;
	set<node *> Lp_zero;

	// Finds Λp
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

	for (ll i = 0; i < L.size(); i++)
	{
		node *Hi = L[i];

		if (Lp.find(Hi->key) == Lp.end())
		{
			NLp.push_back(Hi);
		}
	}
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
	//TODO: update comment
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

		// Finds xi and yi
		ll ti2 = new_vertex(K, P, C, Q, order, rho, inBase, p);
		ll bi2 = new_vertex(K, P, C, Q, order, rho, inBase, p);

		new_dummy_lines.push_back(pr(ti2, bi2));

		vector<node *> H_ancestors = ancestors(Bl, H);
		H->ordering.push_back(ti2);
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

//TODO: update comment
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
				// if NBa[i] in H\B*
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
				// if NBa[i] in H\B*
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

	//-----step3-----
	for (auto itr_z = Lp_zero.begin(); itr_z != Lp_zero.end(); itr_z++)
	{
		node *Hi = *itr_z;
		
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
	}

	// Renames bi',ti' as bud tip of Hi
	for (ll I = 0; I < Lp_positive.size(); I++)
	{
		node *Hi = Lp_positive[I].first;
		if (Hi->normal)
		{
			ll former_tip = Hi->tip;
			ll former_bud = Hi->bud;
			remove_vertex(Hi->bud, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);
			remove_vertex(Hi->tip, Ba, NBa, Ba_sub, NBa_sub, inBase, order, C);

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
			Hi->tip = new_dummy_lines[I].first;
			Hi->bud = new_dummy_lines[I].second;
			Hi->normal = true;
		}
	}

	//-----step4-----
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

	//cout << "Blossoms_in_MaximalNormalBlossom" << endl;

	for (ll I = Blossoms_in_MaximalNormalBlossom.size() - 1; I >= 0; I--)
	{
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
				// if NBa[i] in Hi\B*
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
				// if NBa[i] in Hi\B*
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
		node *Hi = Blossoms_in_MaximalNormalBlossom[I];

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

	//-----step5-----
	vector<pr> new_dummy_line;
	for (ll I = 0; I < Blossoms_in_MaximalNormalBlossom.size(); I++)
	{
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

	//cout << "-----------Augment end-------" << endl;
}

void Search_in_Blossom(Matrix &C, node *H, Tree &Bl, vector<Field> &p, vector<ll> &Ba, vector<ll> &NBa, vector<ll> &order, vector<ll> &rho, vector<ll> &K, vector<bool> &inBase, Matrix &Q, vector<vector<ll>> &P, ll N)
{
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
		Search(C, inBase, Ba, NBa, H_Ba, H_NBa, Q, K, H_tree, Bl, p, order, rho, P, N);

		// Updates dual variables
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
	}

	Expand(C, H, Ba, NBa, H_Ba, H_NBa, inBase, order, Bl);

	//cout << "-------Search_in_Blossom end----------" << endl;
	//cout << endl;
}

vector<ll> matroid_parity(Matrix &A, Field w[])
{
	// Creates a matrix C with given matrix A
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

	//cout << "Base Greedy end" << endl;
	if (overflow)
	{
		return {};
	}
	if (Ba.size() != A.row)
	{
		// There is no base.
		return Ba;
	}

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

	Matrix C = matrix_inverse(A_b) * A_nb;

	num_type C_mx = C.maximum_number();

	maximum_absolute_num = max(maximum_absolute_num, C_mx);

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

	// cout << "step2 start" << endl; 
	while (SourceLine_vector(A.col, Ba, inBase))
	{
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
				vector<node *> AllBlossoms = all_blossoms(Blossom);

				// Executes from the largest blossom to the smallest
				for (ll i = 0; i < AllBlossoms.size(); i++)
				{
					node *Hi = AllBlossoms[i];
					if (Hi->q.isZero() && isMaximal(Blossom, Hi))
					{
						Expand(C, Hi, Ba, NBa, Ba_sub, NBa_sub, inBase, order, Blossom);
						if (overflow)
						{
							return {};
						}
					}
				}
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
				incidence_matrix.X[A[j][l]][edge_num + cnt] = RationalNumber(1, 1);
				incidence_matrix.X[vertex_num + i][edge_num + cnt] = -RationalNumber(-1, 1);
				cnt++;
				w.push_back(RationalNumber(0, 1));
			}
		}
		A[A.size() - 1].push_back(vertex_num + i);
	}

	incidence_matrix.X[vertex_num + additional_terminal_num][edge_num + cnt + 1] = RationalNumber(1, 1);
	incidence_matrix.X[vertex_num + additional_terminal_num + 1][edge_num + cnt + 1] = RationalNumber(-1, 1);

	A.resize(A.size() + 1);
	A[A.size() - 1].push_back(vertex_num + additional_terminal_num);
	w.push_back(RationalNumber(0, 1));
	A.resize(A.size() + 1);
	A[A.size() - 1].push_back(vertex_num + additional_terminal_num + 1);
	w.push_back(RationalNumber(0, 1));

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
}

int main()
{

	int mode;
	cout << "ランダムなら０，接続行列なら1，接続行列の変種なら２を入力" << endl;
	cin >> mode;

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
						filepath = "./test_cases/data_";
						filepath += to_string(weight) + "_" + to_string(element) + "_" + nonzero_prob;
						filepath += "/data_" + to_string(rowsize) + "_" + to_string(rowsize * 2) + "_" + to_string(weight) + "_" + to_string(element) + "_" + nonzero_prob + "_" + to_string(file_num) + ".txt";
						cout << filepath << endl;
					}
					if (mode == 1)
					{
						filepath = "./test_cases/incidence_un_3/incidence_matrix_";
						filepath += to_string(rowsize) + "_" + to_string(rowsize * 3) + "_" + to_string(file_num) + ".txt";
						cout << filepath << endl;
					}
					if (mode == 2)
					{
						filepath = "./test_cases/spath_";
						filepath += to_string(element) + "_0.5/spath_";
						filepath += to_string(element) + "_" + to_string(rowsize) + "_" + to_string(rowsize * 3) + "_" + to_string(file_num) + "_" + "0.5" + ".txt";
					}
					if (mode == 4)
					{
						filepath = "./test_cases/yuukou_inc/yuukou_incidence_matrix_";
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
					Matrix A(stoll(x), stoll(y));
					A.matrix_setZero();
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
						}
						cnt++;
						if (cnt == A.row - 1)
						{
							for (ll i = 0; i < A.col / 2; i++)
							{
								string str;
								ifs >> str;
								ll num = stoll(str);
								w[i].input_ll(num);
							}
							break;
						}
					}

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
