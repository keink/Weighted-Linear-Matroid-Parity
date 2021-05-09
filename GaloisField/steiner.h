#include<iostream>
#include<vector>
#include<algorithm>
#include"Matrix_ov.h"
#include<fstream>
#include <string>

using namespace std;

typedef long long ll;


bool isNum(char a)
{
    if (0 <= a - '0' && a - '0' <= 9)
    {
        return true;
    }
    return false;
}

void input_for_steiner(string filename,vector<ll>& terminals,vector<vector<pair<ll,ll>>>& G)
{

    ifstream ifs;
    ifs.open(filename);
    if (ifs.fail())
    {
        cout << "Fail" << endl;
        return;
    }

    int node_num;
    int edge_num;
    int terminals_num;

    while (!ifs.eof())
    {
        string str;
        ifs >> str;

        if (str == "Nodes")
        {
            string str1;
            ifs >> str1;
            node_num = stoi(str1);
            G.resize(node_num);
        }
        else if (str == "Edges")
        {
            string str1;
            ifs >> str1;
            edge_num = stoi(str1);

            string str2, str3, str4;
            for (int i = 0; i < edge_num; i++)
            {
                ifs >> str1 >> str2 >> str3 >> str4;
                int v1 = stoi(str2);
                int v2 = stoi(str3);
                v1--; v2--;
                int w = stoi(str4);
                G[v1].push_back(make_pair(v2, w));
                G[v2].push_back(make_pair(v1, w));
            }
        }
        else if (str == "Terminals")
        {
            string str1;
            ifs >> str1;
            ifs>>str1;
            terminals_num = stoi(str1);
            string str2;
            for (int i = 0; i < terminals_num; i++)
            {
                ifs >> str1 >> str2;
                int terminal=stoi(str2);
                terminal--;
                terminals.push_back(terminal);
            }
        }
    }
    return;
}

void distance_between_allpairs(vector<vector<pair<ll,ll>>>& G, vector<vector<ll>>& dist)
{
    int n=G.size();

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
        {
            if(i==j){dist[i][j]=0;}
            else{dist[i][j]=INF;}
        }
    }

    
    for(int i=0;i<G.size();i++){
        for(int j=0;j<G[i].size();j++)
        {
            dist[i][G[i][j].first]=G[i][j].second;
        }
    }
    
    
    for(int k=0;k<n;k++)
    {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                dist[i][j]=min(dist[i][j],dist[i][k]+dist[k][j]);
            }
        }
    }
}

//↓dreyfus and wagner
//algothm of dreyfus and wagner is for terminal sets with arbitrary terminals,
//but below is only for two or three terminals
ll minimum_steiner_tree(vector<long long>& terminals,vector<vector<pair<ll,ll>>>& G,vector<vector<long long>>& dist)
{
    int n=G.size();
 
    if(terminals.size()==2)
    {
        return dist[terminals[0]][terminals[1]];
    }
    else{
        ll ret=INF;
        for(int i=0;i<n;i++)
        {
            ret=min(ret,dist[terminals[0]][i]+dist[terminals[1]][i]+dist[terminals[2]][i]);
        }
        return ret;
    }
}

pair<Matrix,pair<vector<ll>,vector<ll>>> create_matrix_for_steiner(vector<long long>& terminals,vector<vector<pair<ll,ll>>>& G)
{
    vector<vector<ll>> dist;
    dist.resize(G.size());
    for(int i=0;i<G.size();i++)
    {
        dist[i].resize(G.size());
    }

    distance_between_allpairs(G,dist);
    vector<pair<long long, long long> > edgeWithTwoVertices;
    vector<pair<long long, pair<long long, long long> > > edgeWithThreeVertices;

    vector<ll> w_,w;
    int ts=terminals.size();
    int w_size(ts*(ts-1)/2);
    int wsize=(ts*(ts-1)*(ts-2)/6);

    Matrix A(terminals.size(),wsize*2+w_size);
    A.matrix_setZero();

    int cnt=0;
    for(int i=0;i<terminals.size();i++){
        for(int j=i+1;j<terminals.size();j++)
        {
            edgeWithTwoVertices.push_back(make_pair(i,j));
            vector<ll> v={i,j};
            w_.push_back(minimum_steiner_tree(v,G,dist));
            A.X[i][wsize*2+cnt].setOne();
            A.X[j][wsize*2+cnt].setOne();
            cnt++;
        }
    }
    cnt=0;
    for(int i=0;i<terminals.size();i++){
        for(int j=i+1;j<terminals.size();j++){
            for(int k=j+1;k<terminals.size();k++){
                edgeWithThreeVertices.push_back(make_pair(i,make_pair(j,k)));
                vector<ll> v={i,j,k};
                w.push_back(minimum_steiner_tree(v,G,dist));
                A.X[i][cnt].setOne();
                A.X[j][cnt].setOne();
                A.X[i][cnt+1].setOne();
                A.X[k][cnt+1].setOne();
                cnt+=2;
            }
        }
    }
    
    return make_pair(A,make_pair(w,w_));
}