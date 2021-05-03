#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <queue>

typedef long long ll;

class node
{
public:
	class node *parent;
	class node *fchild; // First child
	class node *next;
	ll key;
	ll label;
	Field q;
	bool normal;
	std::vector<std::vector<ll>> routing;
	std::vector<ll> ordering;
	ll tip;
	ll bud;

	node()
	{
		fchild = NULL;
		next = NULL;
		key = 0;
		label = 0;
		q.setZero();
		tip = 0;
		bud = 0;
	}
	~node()
	{
		for (ll i = 0; i < routing.size(); i++)
		{
			routing[i].clear();
			routing[i].shrink_to_fit();
		}
		routing.clear();
		routing.shrink_to_fit();

		ordering.clear();
		ordering.shrink_to_fit();
	}
};

typedef std::pair<node *, node *> p1;
typedef std::pair<node *, p1> pp;

class Tree
{
public:
	node *root;
	ll num;

	Tree()
	{
		node *r = new node();
		num = 0;
		r->key = -1;
		r->next = NULL;
		r->parent = NULL;
		root = r;
	}

	// root is a dummy node
	Tree(ll n)
	{
		node *r = new node();
		num = n;

		r->key = -1;
		r->next = NULL;
		r->parent = NULL;
		root = r;

		node *x = new node();
		x->key = 0;
		r->fchild = x;
		x->parent = r;
		x->ordering.push_back(x->key);
		for (ll i = 1; i < n; i++)
		{
			node *y = new node();
			y->key = i;
			y->ordering.push_back(i);
			x->next = y;
			y->parent = r;
			x = y;
		}
		x->next = NULL;
	}

	~Tree()
	{
	}
};

bool isLeaf(Tree &T, node *x)
{
	if (x->fchild == NULL)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool isMaximal(Tree &T, node *x)
{
	if (x->parent == T.root)
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<node *> maximal_blossom(Tree &T)
{
	node *x;
	std::vector<node *> v;

	x = (T.root)->fchild;
	while (x != NULL)
	{
		if (!isLeaf(T, x))
		{
			v.push_back(x);
		}
		x = x->next;
	}
	return v;
}

std::vector<node *> descendant(Tree &T, node *H)
{
	node *x = H;
	std::queue<node *> que;
	std::vector<node *> ret;
	que.push(x);

	while (!que.empty())
	{
		node *x = que.front();
		node *y = x->fchild;
		que.pop();
		while (y != NULL)
		{
			que.push(y);
			if (!isLeaf(T, y))
			{
				ret.push_back(y);
			}
			y = y->next;
		}
	}
	return ret;
}

std::vector<node *> all_blossoms(Tree &T)
{
	node *x = T.root;

	std::queue<node *> que;
	std::vector<node *> ret;
	que.push(x);

	while (!que.empty())
	{
		node *x = que.front();
		node *y = x->fchild;
		que.pop();
		while (y != NULL)
		{
			que.push(y);
			if (!isLeaf(T, y))
			{
				ret.push_back(y);
			}
			y = y->next;
		}
	}
	return ret;
}

std::vector<node *> blossoms_in_MaximalNormalBlossom(Tree &T)
{
	node *x = T.root;

	std::queue<node *> que;
	std::vector<node *> ret;

	que.push(T.root);

	while (!que.empty())
	{
		node *x = que.front();
		node *y = x->fchild;
		que.pop();
		while (y != NULL)
		{
			if (y->normal)
			{
				que.push(y);
				if (!isLeaf(T, y))
				{
					ret.push_back(y);
				}
			}
			y = y->next;
		}
	}
	return ret;
}

std::vector<node *> ancestors(Tree &T, node *x)
{

	std::vector<node *> ret;
	x = x->parent;
	while (x != T.root)
	{
		ret.push_back(x);

		x = x->parent;
	}
	return ret;
}

pp tree_search(Tree &T, ll k)
{
	std::queue<node *> que;
	que.push(T.root);

	while (!que.empty())
	{
		node *x = que.front();
		node *y = x->fchild;
		node *prev = NULL;
		que.pop();
		while (y != NULL)
		{
			if (y->key == k)
			{
				return pp(y, p1(x, prev));
			}
			prev = y;
			que.push(y);
			y = y->next;
		}
	}
	return pp(NULL, p1(NULL, NULL));
}

void blossom_insert(Tree &T, ll N, std::set<ll> s, Field qh)
{
	pp x = tree_search(T, *(s.begin()));

	node *p = x.second.first;
	node *q = p->fchild;

	// Inserts y
	node *y = new node();
	y->key = N + T.num;
	y->q = qh;
	y->parent = p;

	p->fchild = y;

	// Divides children into two groups: 
	// those that will be child of y and those that will be a sibling of y
	ll cnt = 0;
	node *prev1 = new node();
	node *prev2 = y;

	while (q != NULL)
	{
		if (s.find(q->key) != s.end())
		{
			// Makes q a child
			if (cnt == 0)
			{
				q->parent = y;
				y->fchild = q;
				prev1 = q;
			}
			else
			{
				q->parent = y;
				prev1->next = q;
				prev1 = q;
			}
			cnt++;
		}
		else
		{
			// Makes q a sibling
			prev2->next = q;
			prev2 = q;
		}
		q = q->next;
	}
	prev1->next = NULL;
	prev2->next = NULL;
}

void tree_delete(Tree &T, ll k)
{
	if (k == -1)
	{
		return;
	}

	pp x = tree_search(T, k);
	node *parent_node = x.second.first;
	node *first_child_node = x.first->fchild;
	node *next_node = x.first->next;
	node *previous_node = x.second.second;

	if (previous_node != NULL)
	{
		previous_node->next = first_child_node;
		node *y = first_child_node;
		node *z = previous_node;
		while (y != NULL)
		{
			y->parent = parent_node;
			z = y;
			y = y->next;
		}
		if (z != NULL)
		{
			z->next = next_node;
		}
	}
	else
	{
		if (first_child_node == NULL)
		{
			parent_node->fchild = next_node;
		}
		else
		{
			parent_node->fchild = first_child_node;
			node *y = first_child_node;
			node *z;
			while (y != NULL)
			{
				y->parent = parent_node;
				z = y;
				y = y->next;
			}
			z->next = next_node;
		}
	}
	delete x.first;
}

void blossom_initialize(Tree &T, std::vector<ll> &K)
{
	std::queue<node *> que;
	que.push(T.root);

	while (!que.empty())
	{
		node *x = que.front();
		node *y = x->fchild;
		que.pop();
		while (y != NULL)
		{

			if (x == T.root && !isLeaf(T, y))
			{
				for (ll i = 0; i != y->ordering.size(); i++)
				{
					K[y->ordering[i]] = y->key;
				}
				if (y->normal)
				{
					K[y->bud] = y->key;
				}
			}

			que.push(y);
			y->label = 0;
			y = y->next;
		}
	}
}
