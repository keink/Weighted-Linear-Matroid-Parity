# Description
This is an implementation of an algorithm for weighted linear matroid parity problem, proposed by Iwata and Kobayashi [1].

# Requirements
[Boost Multiprecision Library](https://www.boost.org/doc/libs/1_75_0/libs/multiprecision/doc/html/index.html)

# Usage
You can switch types of matrices by standard input.
Mode 0, 1, 2, 4 respectively correspond to integer random matrices, incidence matrices of undirected graphs, matrices for S-path problem and incidence matrices of directed graphs.

## rationalNumberField/WeightedLinearMatroidParity.cpp 
the algorithm over rational number field

## GaloisField/WeightedLinearMatroidParity.cpp
the algorithm over Galois Field
Also it supports partial pairng (Parts of columns are paired).

## Test cases
Each test case starts with the number of rows and the number of columns, and elements of matrix and weights of lines follow.
For instance, if you want to run pro the following matrix and weights of lines are (2, 3),
https://latex.codecogs.com/gif.latex?\begin{pmatrix}&space;0&space;&&space;1&space;&&space;2&space;&&space;3&space;\\&space;-1&space;&&space;0&space;&&space;1&space;&&space;2&space;\\&space;-2&space;&&space;-1&space;&&space;0&space;&&space;1&space;\\&space;-3&space;&&space;-2&space;&&space;-1&space;&&space;0&space;\end{pmatrix}
then test case should be written as below.
***
4 4
0 1 2 3
-1 0 1 2
-2 -1 0 1
-3 -2 -1 0
1 2
***

## Output
For each test case, programs output the following terms:  <br>
- initial base found by the algorithm<br>
- how many times each procedure of the algorithm executed,<br>
- maximum absolute number appeared in the target matrix while executing the algorithm, <br>
- solution, <br>
- optimal value,<br>
- running time.

# References
[1] S. Iwata and Y. Kobayashi: A weighted linear matroid parity algorithm, Proceedings of the 49th ACM Symposium on Theory of Computing (STOC), 2018, pp. 264-276

