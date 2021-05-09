#pragma once
#include <iostream>
#include <vector>
#include "Field_gf.h"
#include "RationalNumber_gf.h"
#define Field RationalNumber
#include "Matrix2_gf.h"


class Matrix : public Field
{
  public:
	ll row;
	ll col;
	std::vector<std::vector<Field>> X;

	Matrix(ll r, ll c)
	{
		row = r;
		col = c;
		X.resize(row);
		for (ll i = 0; i < row; i++)
		{
			X[i].resize(col);
		}

		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				if (i == j)
					X[i][j].setOne();
				else
					X[i][j].setZero();
			}
		}
	}

	Matrix()
	{
		row = 0;
		col = 0;
	}

	~Matrix()
	{
		for (ll i = 0; i < row; i++)
		{
			X[i].clear();
			X[i].shrink_to_fit();
		}
		X.clear();
		X.shrink_to_fit();
	}

	Matrix operator -() const{
		Matrix A(row,col);
		for(ll i=0;i<row;i++){
			for(ll j=0;j<col;j++){
				A.X[i][j]=-X[i][j];
			}
		}
		return A;
	}

	void input_matrix()
	{
		std::cin >> row >> col;

		X.resize(row);
		for (ll i = 0; i < row; i++)
		{
			X[i].resize(col);
		}

		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				X[i][j].input();
			}
		}
	}

	void input_matrix_int()
	{
		std::cin>>row>>col;
		X.resize(row);
		for (ll i = 0; i < row; i++)
		{
			X[i].resize(col);
		}
		ll x;
		for(ll i=0;i<row;i++){
			for(ll j=0;j<col;j++){
				std::cin>>x;
				X[i][j]=RationalNumber(x,1);
			}
		}

	}

	void output_matrix()
	{
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				X[i][j].output();
				std::cout << " ";
			}
			std::cout << std::endl;
		}
	}

	void output_matrix_num()
	{
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				X[i][j].output_num();
				std::cout << " ";
			}
			std::cout << std::endl;
		}
	}

	void insert()
	{
	}

	void row_swap(ll i, ll j)
	{
		if (i != j)
		{
			for (ll k = 0; k < col; k++)
			{
				std::swap(X[i][k], X[j][k]);
			}
		}
	}

	void column_swap(ll i, ll j)
	{
		if (i != j)
		{
			for (ll k = 0; k < row; k++)
			{
				std::swap(X[k][i], X[k][j]);
			}
		}
	}

	void row_multiply(ll i, Field c)
	{
		for (ll k = 0; k < col; k++)
		{
			X[i][k] *= c;
		}
	}

	void column_multiply(ll j, Field c)
	{
		for (ll k = 0; k < row; k++)
		{
			X[k][j] *= c;
		}
	}

	void row_plus(ll i, ll j, Field c)
	{
		for (ll k = 0; k < col; k++)
		{
			X[i][k] += c * X[j][k];
		}
	}

	void column_plus(ll i, ll j, Field c)
	{
		for (ll k = 0; k < row; k++)
		{
			X[k][i] += c * X[k][j];
		}
	}

	Matrix operator*(const Matrix &A)
	{
		Matrix B = Matrix(row, A.col);

		for(ll i=0;i<B.row;i++){
			for(ll j=0;j<B.col;j++){
				B.X[i][j].setZero();
			}
		}
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < A.col; j++)
			{
				for (ll k = 0; k < col; k++)
				{
					B.X[i][j] += X[i][k] * A.X[k][j];
				}
			}
		}
		return B;
	}
	Matrix operator-(const Matrix& A)
	{
		Matrix B=Matrix(row,col);
		for(ll i=0;i<row;i++){
			for(ll j=0;j<col;j++){
				B.X[i][j]=X[i][j]-A.X[i][j];
			}
		}
		return B;
	}

	void matrix_setZero(){
		for(ll i=0;i<row;i++){
			for(ll j=0;j<col;j++){
				X[i][j].setZero();
			}
		}
	}

	Matrix sub_matrix(std::vector<ll> A, std::vector<ll> B)
	{
		Matrix S(A.size(), B.size());
		for (ll i = 0; i < A.size(); i++)
		{
			for (ll j = 0; j < B.size(); j++)
			{
				S.X[i][j] = X[A[i]][B[j]];
			}
		}
		return S;
	}

	num_type maximum_number()
    {
        num_type gamma=0;
        for(ll i=0;i<row;i++){
            for(ll j=0;j<col;j++){
                gamma=std::max(gamma,X[i][j].max_num());
            }
        }
        return gamma;
    }

	Matrix2 RNtoGF()
	{
		Matrix2 B(row,col);
		for(ll i=0;i<row;i++){
			for(ll j=0;j<col;j++)
			{
				B.X[i][j]=GF(X[i][j].numerator());
			}
		}

	return B;
	}

};

Matrix matrix_inverse(const Matrix& A)
{
	Matrix B=A;
	ll c = A.col;
	ll r = A.row;

	if (r != c)
	{
		std::cout << "The size of this matrix is " << A.row << "*" << A.col << ", so we cannot compute inverse matrix." << std::endl;
	}

	ll n = c;

	Matrix A_inverse(n, n);

	for (ll j = 0; j < n; j++)
	{
		Field mx;
		mx.setZero();
		ll pivot = j;
		for (ll i = j; i < n; i++)
		{
			if (mx < B.X[i][j] || mx < -B.X[i][j])
			{
				mx = std::max(B.X[i][j], -B.X[i][j]);
				pivot = i;
			}
		}
		if (mx.isZero())
		{
			std::cout<<"A is not nonsingular"<<std::endl;
			exit(EXIT_FAILURE);
		}

		B.row_swap(j, pivot);
		A_inverse.row_swap(j, pivot);

		Field a = B.X[j][j].inverse();

		// Multiplies by a so that A[j][j]=1
		B.row_multiply(j, a);
		A_inverse.row_multiply(j, a);

		for (ll i = 0; i < n; i++)
		{
			if (i != j)
			{
				Field b = B.X[i][j];
				B.row_plus(i, j, -b);
				A_inverse.row_plus(i, j, -b);
			}
		}
	}
	return A_inverse;
}

Matrix toINT(Matrix &A)
{
	Matrix B(A.row,A.col);

	num_type L=1;
	for(ll i=0;i<A.row;i++){
		for(ll j=0;j<A.col;j++){
			L=lcm(L,A.X[i][j].denominator());
		}
	}

	for(ll i=0;i<A.row;i++){
		for(ll j=0;j<A.col;j++){
			B.X[i][j]=A.X[i][j]*RationalNumber(L,1);
		}
	}
	return B;
}