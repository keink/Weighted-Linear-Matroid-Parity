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
		//std::cout<<"overflow"<<std::endl;
		//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
		//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
		overflow = true;
		//exit(EXIT_FAILURE);
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
		//std::cout<<"overflow"<<std::endl;
		//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
		//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
		overflow = true;
		//exit(EXIT_FAILURE);
	}
	ll k = c * d;
	if (overflow_check_product(k, g))
	{
		//std::cout<<"overflow"<<std::endl;
		//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
		//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
		overflow = true;
		//exit(EXIT_FAILURE);
	}
	return k * g;
}

class RationalNumber
{
	num_type num; //bunshi
	num_type den; //bunbo

public:
	RationalNumber(num_type numerator, num_type denominator)
	{
		/* if(numerator>LLONG_MAX||denominator>LLONG_MAX){
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);
		} */
		if (denominator == 0)
		{
			//std::cout<<"num is "<<numerator<<" den is "<<denominator<<std::endl;
			//std::cout << "Error: denominator is 0" << std::endl;
			/* if(overflow){
				//std::cout<<"!!!!!!!!!!!!!!!!!"<<std::endl;
				//std::cout<<"overflow"<<std::endl;
				denominator=1;
				numerator=1;
				overflow=true;
			}
			else{ */
			std::cout << "Error: denominator is 0" << std::endl;
			exit(EXIT_FAILURE);
			//}
		}
		if (denominator < 0)
		{
			/* if(numerator==LLONG_MIN||denominator==LLONG_MIN){
				//std::cout<<"overflow"<<std::endl;
				//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
				//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
				overflow=true;
				//exit(EXIT_FAILURE);
			} */
			numerator = -numerator;
			denominator = -denominator;
		}
		num_type g = boost::multiprecision::gcd(numerator, denominator);

		num = numerator / g;
		den = denominator / g;
	}

	/* 	RationalNumber(ll numerator){
		num=numerator; den=1;
	} */

	RationalNumber()
	{
		num = 0;
		den = 1;
	}

	~RationalNumber()
	{
		//STD::cout<<col<<" "<<row<<std::endl;
		//std::cout<<"Rational Number デストラクタ"<<std::endl;
	}

	RationalNumber operator+() const { return *this; }
	RationalNumber operator-() const
	{
		/* if(num==LLONG_MIN){
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);
		} */
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
		/* if(overflow_check_product(num,lc/den)||overflow_check_product(x.num,lc/x.den))
		{
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);
			return true;
		} */
		return (num * (lc / den)) < (x.num * (lc / x.den));
	}

	RationalNumber &operator=(const RationalNumber &x)
	{
		//std::cout<<den<<std::endl;
		den = x.den;
		num = x.num;
		//std::cout<<"c"<<std::endl;
		return *this;
	}

	RationalNumber operator+(const RationalNumber &x) const
	{
		//ll de = den * x.den;
		//ll nu = num * x.den + den * x.num;

		num_type de = boost::multiprecision::lcm(den, x.den);
		/* if(de==LLONG_MIN&&den==-1){
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);
		}
		if(de==LLONG_MIN&&x.den==-1){
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);	
		}
		if(overflow_check_product(de/den,num)||overflow_check_product(de/x.den,x.num))
		{
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;
			//exit(EXIT_FAILURE);
		}
		if(overflow_check_plus((de/den)*num,(de/x.den)*x.num))
		{
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;	
			overflow=true;//exit(EXIT_FAILURE);
		} */
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

		/* if(overflow_check_product(nu1,nu2)||overflow_check_product(de1,de2)){
			//std::cout<<"overflow"<<std::endl;
			//std::cout<<"the number of times of excuting Search/Blossom/Graft/DualUpdate/Augment"<<std::endl;
			//std::cout<<Search_num<<" "<<Blossom_num<<" "<<Graft_num<<" "<<DualUpdate_num<<" "<<Augment_num<<std::endl;
			overflow=true;//exit(EXIT_FAILURE);
		} */

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

	/*
	inline std::ostream& operator<<(std::ostream& s)
	{
		s<<num<<"/"<<den;
		return s;
	}
	////////
	inline std::istream& operator>>(std::istream& s)
	{
		s>>den;
		return s;
	}
	*/

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
		/*
		if(num%2==0)
		{
			a.num/=2;
		}
		else{
			a.den*=2;
		}
		*/
		return a * RationalNumber(1, 2);
	}

	RationalNumber Double() const
	{
		RationalNumber a = RationalNumber(num, den);
		/*
		if(den%2==0){
			a.den/=2;
		}
		else{
			a.num*=2;
		}
		*/

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
		//std::cout<<num<<" "<<den<<std::endl;
		RationalNumber x = (*this).abs();

		//std::cout<<"mm"<<std::endl;
		//std::cout<<x.num<<" "<<x.den<<std::endl;
		return std::max(x.num, x.den);
	}
};