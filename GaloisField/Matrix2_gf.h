#pragma once
#include <iostream>
#include <vector>
#include "Field_gf.h"
#include "RationalNumber_gf.h"
#include "GF_gf.h"
#define Field2 GF

class Matrix2 : public Field2
{
  public:
	ll row;
	ll col;
	std::vector<std::vector<Field2>> X;

	Matrix2(ll r, ll c)
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

	Matrix2()
	{
		row = 0;
		col = 0;
	}

	~Matrix2()
	{
		for (ll i = 0; i < row; i++)
		{
			X[i].clear();
			X[i].shrink_to_fit();
		}
		X.clear();
		X.shrink_to_fit();
	}

	Matrix2 operator -() const{
		Matrix2 A(row,col);
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

	void row_multiply(ll i, Field2 c)
	{
		for (ll k = 0; k < col; k++)
		{
			X[i][k] *= c;
		}
	}

	void column_multiply(ll j, Field2 c)
	{
		for (ll k = 0; k < row; k++)
		{
			X[k][j] *= c;
		}
	}

	void row_plus(ll i, ll j, Field2 c)
	{
		for (ll k = 0; k < col; k++)
		{
			X[i][k] += c * X[j][k];
		}
	}

	void column_plus(ll i, ll j, Field2 c)
	{
		for (ll k = 0; k < row; k++)
		{
			X[k][i] += c * X[k][j];
		}
	}

	Matrix2 operator*(const Matrix2 &A)
	{
		Matrix2 B = Matrix2(row, A.col);

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
	Matrix2 operator-(const Matrix2& A)
	{
		Matrix2 B=Matrix2(row,col);
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

	Matrix2 sub_matrix(std::vector<ll> A, std::vector<ll> B)
	{
		Matrix2 S(A.size(), B.size());
		for (ll i = 0; i < A.size(); i++)
		{
			for (ll j = 0; j < B.size(); j++)
			{
				S.X[i][j] = X[A[i]][B[j]];
			}
		}
		return S;
	}
	
};

Matrix2 matrix2_inverse(const Matrix2& A)
{
	Matrix2 B=A;
	ll c = A.col;
	ll r = A.row;

	if (r != c)
	{
		std::cout << "The size of this matrix is " << A.row << "*" << A.col << ", so we cannot compute an inverse matrix." << std::endl;
	}

	ll n = c;

	Matrix2 A_inverse(n, n);

	for (ll j = 0; j < n; j++)
	{
		Field2 mx;
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

		Field2 a = B.X[j][j].inverse();

		// Multiplies by a so that A[j][j]=1
		B.row_multiply(j, a);
		A_inverse.row_multiply(j, a);

		for (ll i = 0; i < n; i++)
		{
			if (i != j)
			{
				Field2 b = B.X[i][j];
				B.row_plus(i, j, -b);
				A_inverse.row_plus(i, j, -b);
			}
		}
	}
	return A_inverse;
}