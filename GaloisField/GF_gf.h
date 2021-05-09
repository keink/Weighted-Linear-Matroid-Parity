#pragma once
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <numeric>
#include "Field_gf.h"
#include <boost/multiprecision/cpp_int.hpp>

typedef long long ll;

extern ll prime_number;

ll mod(ll x, ll p){
    if(x>=0){
        return x%p;
    }
    else{
        return (x+((-x)/p+1)*p)%p;
    }
}

ll extentedEuglid(ll a,ll b, ll& x, ll& y)
{
    if(b!=0){
        ll d=extentedEuglid(b,a%b,y,x);
        y-=(a/b)*x;
        return d;
    }
    
        x=1;y=0;

    return a;
}

class GF: public Field
{
    
public:
ll num;
    GF(ll n)
    {
        num=mod(n,prime_number);
    }

    GF()
    {
        num=0;
    }

    ~GF()
    {

    }

    GF operator+()const{return *this;}
    GF operator-()const{return GF(prime_number-num);}

    ///////
    GF inverse()const{
        ll x; ll y;
        ll d=extentedEuglid(num,prime_number,x,y);
        return GF(x);
    }

    bool isZero()
    {
        return mod(num,prime_number)==0;
    }

    bool operator==(const GF &x)
    {
        return mod(num,prime_number)==mod(x.num,prime_number);
    }

    bool operator!=(const GF & x){
        return !(*this==x);
    }

    bool operator<(const GF& x){
        return (num<x.num);
    }

    GF& operator=(const GF& x)
    {
        num=x.num;
        return *this;
    }

    GF operator+(const GF& x)
    {
        return GF(num+x.num);
    }

    GF &operator+=(const GF& x){
        *this=*this+x;
    }

    GF operator-(const GF& x){
        return *this+(-x);
    }

    GF &operator-=(const GF& x)
    {
        *this=*this-x;
    }

    GF operator*(const GF& x)
    {
        return GF(num*x.num);
    }

    GF &operator*=(const GF& x)
    {
        *this=*this*x;
    }

    //////////
    GF operator/(const GF& x)
    {
        return (*this)*x.inverse();
    }

    GF &operator/=(const GF& x){
        *this=*this/x;
    }

    void output()
    {
        std::cout<<num;
    }

    void input()
    {
        std::cin>>num;
        num=mod(num,prime_number);
    }

    void setZero()
    {
        num=0;
    }

    void setOne()
    {
        num=1;
    }

    bool operator<(const GF &x) const
    {
        return num<x.num;
    }
};

// Returns a list of all prime numbers which are less than or equal to n
std::vector<ll> eratosthenes(ll n){
    std::vector<bool> is_prime;
    std::vector<ll> prime;

    ll N=n*n;
    is_prime.resize(N+1,true);

    is_prime[0]=is_prime[1]=false;

    for(ll i=2;i<=N;i++){
        if(is_prime[i]==true){
            prime.push_back(i);
            for(ll j=2*i;j<=N;j+=i){
                is_prime[j]=false;
            }
        }
        if(prime.size()==n){
            break;
        }
    }
    return prime;

}

// Returns a list of all prime numbers which are less than or equal to n
std::vector<ll> eratosthenes_(ll n){
    std::vector<bool> is_prime;
    std::vector<ll> prime;

    is_prime.resize(n+1,true);

    is_prime[0]=is_prime[1]=false;

    for(ll i=2;i<=n;i++){
        if(is_prime[i]==true){
            prime.push_back(i);
            for(ll j=2*i;j<=n;j+=i){
                is_prime[j]=false;
            }
        }
    }
    return prime;
}