#pragma once
#include <iostream>
#include <vector>
#include "RationalNumber.h"

using Field = RationalNumber;

class Matrix
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
				//std::cout<<i<<"  "<<j<<std::endl;
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
		//std::cout<<row<<" "<<col<<std::endl;
		//std::cout<<"Matrix デストラクタ"<<std::endl;
		for (ll i = 0; i < row; i++)
		{
			X[i].clear();
			X[i].shrink_to_fit();
		}
		X.clear();
		X.shrink_to_fit();
	}

	Matrix operator-() const
	{
		Matrix A(row, col);
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				A.X[i][j] = -X[i][j];
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
				//std::cout<<i<<"  "<<j<<std::endl;
			}
		}
	}

	void input_matrix_int()
	{
		std::cin >> row >> col;
		X.resize(row);
		for (ll i = 0; i < row; i++)
		{
			X[i].resize(col);
		}
		ll x;
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				std::cin >> x;
				X[i][j] = RationalNumber(x, 1);
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

		for (ll i = 0; i < B.row; i++)
		{
			for (ll j = 0; j < B.col; j++)
			{
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
	Matrix operator-(const Matrix &A)
	{
		Matrix B = Matrix(row, col);
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				B.X[i][j] = X[i][j] - A.X[i][j];
			}
		}
		return B;
	}

	void matrix_setZero()
	{
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
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
				//std::cout << i << " " << j << " " << A[i] << " " << B[j] << std::endl;
				//S.X[i][j].output();
				//std::cout << " ";
				//X[A[i]][B[j]].output();
				//std::cout << std::endl;
			}
		}
		return S;
	}

	num_type maximum_number()
	{
		num_type gamma = 0;
		for (ll i = 0; i < row; i++)
		{
			for (ll j = 0; j < col; j++)
			{
				//std::cout<<i<<" "<<j<<std::endl;
				gamma = std::max(gamma, X[i][j].max_num());
				//std::cout<<i<<" "<<j<<std::endl;
			}
		}
		return gamma;
	}
};

Matrix matrix_inverse(const Matrix &A)
{
	Matrix B = A;
	ll c = A.col;
	ll r = A.row;

	if (r != c)
	{
		std::cout << "この行列のサイズは" << A.row << "*" << A.col << "なので逆行列を求められません" << std::endl;

		//exit(EXIT_FAILURE);
	}

	ll n = c;

	Matrix A_inverse(n, n);

	//A_inverse.output_matrix();

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
			std::cout << "A is not nonsingular" << std::endl;
			//A_inverse.matrix_setZero();
			//return A_inverse;
			exit(EXIT_FAILURE);
		}

		B.row_swap(j, pivot);
		A_inverse.row_swap(j, pivot);

		/* 			std::cout<<"aaa"<<std::endl;
			output_matrix();
			std::cout<<"aaa"<<std::endl;
			A_inverse.output_matrix();
			std::cout<<std::endl; */

		Field a = B.X[j][j].inverse();

		//multiply by a so that A[j][j]=1
		B.row_multiply(j, a);
		A_inverse.row_multiply(j, a);

		for (ll i = 0; i < n; i++)
		{
			if (i != j)
			{
				Field b = B.X[i][j];
				B.row_plus(i, j, -b);
				A_inverse.row_plus(i, j, -b);
				/*				for(ll k=0;k<n;k++){
					A.row_plus(k,j,-b);
					A.X[i][k]-=A.X[j][k]*b;
					A_inverse.X[i][k]-=A_inverse.X[j][k]*b;
				}*/
			}
		}
	}
	return A_inverse;
}

/* Matrix GaussJordan(Matrix &A)
{
	ll c = A.col;
	ll r = A.row;

	if (r != c)
		exit(EXIT_FAILURE);

	ll n = c;

	Matrix A_inverse(n, n);
	//A_inverse.output_matrix();

	for (ll j = 0; j < n; j++)
	{
		Field mx;
		mx.setZero();
		ll pivot = j;
		for (ll i = j; i < n; i++)
		{
			if (mx < A.X[i][j] || mx < -A.X[i][j])
			{
				mx = std::max(A.X[i][j], -A.X[i][j]);
				pivot = i;
			}
		}
		if (mx.isZero())
			exit(EXIT_FAILURE);
		A.row_swap(j, pivot);
		A_inverse.row_swap(j, pivot);

		Field a = A.X[j][j].inverse();

		//multiply by a so that A[j][j]=1
		A.row_multiply(j, a);
		A_inverse.row_multiply(j, a);

		for (ll i = 0; i < n; i++)
		{
			if (i != j)
			{
				Field b = A.X[i][j];
				A.row_plus(i, j, -b);
				A_inverse.row_plus(i, j, -b);
				/*				for(ll k=0;k<n;k++){
					A.row_plus(k,j,-b);
					A.X[i][k]-=A.X[j][k]*b;
					A_inverse.X[i][k]-=A_inverse.X[j][k]*b;
				}*/
/*
			}
		}
	}
	return A_inverse;
}
 */
Matrix kronecker(Matrix &A, Matrix &B)
{
	ll m1 = A.row;
	ll n1 = A.col;
	ll m2 = B.row;
	ll n2 = B.col;

	Matrix C(m1 * m2, n1 * n2);

	for (int i = 0; i < m1; i++)
	{
		for (int j = 0; j < n1; j++)
		{
			for (int k = 0; k < m2; k++)
			{
				for (int l = 0; l < n2; l++)
				{
					C.X[i * m2 + k][j * n2 + l] = A.X[i][j] * B.X[k][l];
				}
			}
		}
	}

	return C;
}
