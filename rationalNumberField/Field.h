#pragma once
#include <iostream>
#include <math.h>
typedef long long ll;

class Field
{
  public:
	virtual Field operator+(){};
	virtual Field operator-(){};
	virtual Field inverse(){};
	virtual bool isZero(){};
	virtual bool operator==(const Field &x){};
	virtual bool operator!=(const Field &x){};
	virtual bool operator<(const Field &x) const {};
	virtual Field& operator=(const Field &X){};
	virtual Field operator+(const Field &x){};
	virtual Field &operator+=(const Field &x){};
	virtual Field operator-(const Field &x){};
	virtual Field &operator-=(const Field &x){};
	virtual Field operator*(const Field &x){};
	virtual Field &operator*=(const Field &x){};
	virtual Field operator/(const Field &x){};
	virtual Field &operator/=(const Field &x){};
	//virtual std::ostream& operator<<(std::ostream& s);
	//virtual std::istream& operator>>(std::istream& s);
	virtual void output(){};
	virtual void input(){};

	virtual void setZero(){};
	virtual void setOne(){};
	virtual Field divideByTwo(){};
	virtual Field Double(){};
	virtual Field abs(){};

	//virtual Field ll_to_RN(ll x){};
};