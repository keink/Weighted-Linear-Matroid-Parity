#pragma once
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <numeric>
#include <climits>
#include <boost/multiprecision/cpp_int.hpp>

typedef long long ll;
extern ll Blossom_num;
extern ll Graft_num;
extern ll Augment_num;
extern ll Search_num;
extern ll DualUpdate_num;
extern bool overflow;

typedef boost::multiprecision::cpp_int num_type;

bool overflow_check_plus(ll a, ll b)
{
	if (a >= 0 && b >= 0)
	{
		if (a <= LLONG_MAX - b)
		{
			return false;
		}
		return false;
	}
	else if (a < 0 && b >= 0)
	{
		return false;
	}
	if (a < 0 && b < 0)
	{
		if (a >= LLONG_MIN - b)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	return false;
}

bool overflow_check_product(ll a, ll b)
{
	if (a == 0 || b == 0)
	{
		return false;
	}
	if (a > 0 && b > 0)
	{
		if (a <= LLONG_MAX / b)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else if (a < 0 && b < 0)
	{
		if (-a <= LLONG_MAX / (-b))
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else if (a > 0 && b < 0)
	{
		if (b >= LLONG_MIN / a)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	if (a >= LLONG_MIN / b)
	{
		return false;
	}
	return true;
}

ll gcd(ll a, ll b)
{
	if (a == LLONG_MIN || b == LLONG_MIN)
	{
		overflow = true;
	}
	a = abs(a);
	b = abs(b);
	if (b == 0)
	{
		return a;
	}
	else
	{
		return gcd(b, a % b);
	}
}

ll lcm(ll a, ll b)
{
	ll g = gcd(a, b);
	ll c = a / g;
	ll d = b / g;
	if (overflow_check_product(c, d))
	{
		overflow = true;
	}
	ll k = c * d;
	if (overflow_check_product(k, g))
	{
		overflow = true;
	}
	return k * g;
}

class RationalNumber
{
	num_type num; // numerator
	num_type den; // denominator

public:
	RationalNumber(num_type numerator, num_type denominator)
	{
		if (denominator == 0)
		{
			std::cout << "Error: denominator is 0" << std::endl;
			exit(EXIT_FAILURE);
		}
		if (denominator < 0)
		{
			numerator = -numerator;
			denominator = -denominator;
		}
		num_type g = boost::multiprecision::gcd(numerator, denominator);

		num = numerator / g;
		den = denominator / g;
	}

	RationalNumber()
	{
		num = 0;
		den = 1;
	}

	~RationalNumber()
	{
	}

	RationalNumber operator+() const { return *this; }
	RationalNumber operator-() const
	{
		return RationalNumber(-num, den);
	}

	RationalNumber inverse() const
	{
		return RationalNumber(den, num);
	}

	bool isZero()
	{
		return num == 0;
	}

	bool operator==(const RationalNumber &x)
	{
		return num == x.num && den == x.den;
	}

	bool operator!=(const RationalNumber &x)
	{
		return !(*this == x);
	}

	bool operator<(const RationalNumber &x) const
	{
		num_type lc = boost::multiprecision::lcm(den, x.den);
		return (num * (lc / den)) < (x.num * (lc / x.den));
	}

	RationalNumber &operator=(const RationalNumber &x)
	{
		den = x.den;
		num = x.num;
		return *this;
	}

	RationalNumber operator+(const RationalNumber &x) const
	{
		num_type de = boost::multiprecision::lcm(den, x.den);
		num_type nu = (de / den) * num + (de / x.den) * x.num;

		return RationalNumber(nu, de);
	}

	RationalNumber &operator+=(const RationalNumber &x)
	{
		*this = *this + x;
		return *this;
	}

	RationalNumber operator-(const RationalNumber &x) const
	{
		return *this + (-x);
	}

	RationalNumber &operator-=(const RationalNumber &x)
	{
		*this = *this - x;
		return *this;
	}

	RationalNumber operator*(const RationalNumber &x)
	{
		num_type nu1 = num;
		num_type de1 = den;
		num_type nu2 = x.num;
		num_type de2 = x.den;

		num_type g1 = boost::multiprecision::gcd(nu1, de2);

		nu1 = nu1 / g1;
		de2 = de2 / g1;

		num_type g2 = boost::multiprecision::gcd(nu2, de1);
		nu2 = nu2 / g2;
		de1 = de1 / g2;

		return RationalNumber(nu1 * nu2, de1 * de2);
	}

	RationalNumber &operator*=(const RationalNumber &x)
	{
		*this = (*this) * x;
		return *this;
	}

	RationalNumber operator/(const RationalNumber &x)
	{
		return (*this) * (x.inverse());
	}

	RationalNumber &operator/=(const RationalNumber &x)
	{
		*this = *this / x;
		return *this;
	}

	void output()
	{
		std::cout << num << "/" << den;
	}

	void output_num()
	{
		std::cout << num;
	}

	void input()
	{
		std::string s;
		std::cin >> s;
		num_type n = 0;
		ll i = 0;
		while (s[i] != '/')
		{
			n *= 10;
			n += s[i] - '0';
			i++;
		}
		i++;

		num_type d = 0;
		while (i < s.length())
		{
			d *= 10;
			d += s[i] - '0';
			i++;
		}
		num = n;
		den = d;
	}

	void input_ll(ll a)
	{
		num = a;
		den = 1;
	}

	void setZero()
	{
		num = 0;
		den = 1;
	}

	void setOne()
	{
		num = 1;
		den = 1;
	}

	RationalNumber divideByTwo() const
	{
		RationalNumber a = RationalNumber(num, den);
		return a * RationalNumber(1, 2);
	}

	RationalNumber Double() const
	{
		RationalNumber a = RationalNumber(num, den);
		return a * RationalNumber(2, 1);
	}

	RationalNumber abs() const
	{
		if (num >= 0 && den > 0)
		{
			return RationalNumber(num, den);
		}
		else if (num < 0 && den > 0)
		{
			return RationalNumber(-num, den);
		}
		else if (num >= 0 && den < 0)
		{
			return RationalNumber(num, -den);
		}
		return RationalNumber(-num, -den);
	}

	num_type max_num() const
	{
		RationalNumber x = (*this).abs();
		return std::max(x.num, x.den);
	}
};