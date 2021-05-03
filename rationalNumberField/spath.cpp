#include <iostream>
#include <cstdio>
#include <math.h>
#include <string>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <vector>
#include <set>
#include <fstream>
#include <map>
#include <queue>
#include <time.h>
#include <stack>

using namespace std;
typedef pair<int, int> p;
typedef long long ll;

ll blossom_number;
ll vertex_number;
ll ordering_number;

ll Blossom_num = 0;
ll Graft_num = 0;
ll Augment_num = 0;
ll Search_num = 0;
ll DualUpdate_num = 0;

ll maximum_absolute_num = 0;
bool overflow = false;

#include "Matrix_ov.h"

pair<Matrix, vector<RationalNumber>> create_matrix_for_spath(vector<pair<pair<int, int>, RationalNumber>> edges, int k, int vertex_num, int edge_num, vector<vector<int>> A, int terminal_num)
{
    vector<RationalNumber> w;
    int additional_terminal_num = terminal_num - 2 * k;

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

    incidence_matrix.output_matrix();

    incidence_matrix = kronecker(incidence_matrix, E2);

    incidence_matrix.output_matrix();

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
    incidence_matrix.output_matrix();

    return make_pair(incidence_matrix, w);
}

int main()
{
    vector<pair<pair<int, int>, RationalNumber>> edges;
    RationalNumber r = RationalNumber(1, 1);
    edges.push_back(make_pair(make_pair(0, 2), r));
    edges.push_back(make_pair(make_pair(1, 2), r));
    edges.push_back(make_pair(make_pair(0, 3), r));
    edges.push_back(make_pair(make_pair(1, 3), r));
    edges.push_back(make_pair(make_pair(1, 4), r));
    edges.push_back(make_pair(make_pair(2, 3), r));
    edges.push_back(make_pair(make_pair(3, 4), r));
    edges.push_back(make_pair(make_pair(0, 4), r));
    int k = 1;
    int vertex_num = 5;
    int edge_num = 8;
    int terminal_num = 3;
    vector<vector<int>> A;
    A.resize(3);
    A[0].push_back(0);
    A[1].push_back(1);
    A[2].push_back(2);

    pair<Matrix, vector<RationalNumber>> p = create_matrix_for_spath(edges, k, vertex_num, edge_num, A, terminal_num);
    RationalNumber *w = new RationalNumber[p.first.col / 2];
    for (int i = 0; i < p.second.size(); i++)
    {
        w[i] = p.second[i];
    }

    return 0;
}