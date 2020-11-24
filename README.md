# Description
This is an implementation of an algorithm for weighted linear matroid parity problem, proposed by Iwata and Kobayashi [1].

# Requirements
Boost Multiprecision Library

# Usage
You can switch types of matrices by standard input.
Mode 0, 1, 2, 4 respectively correspond to integer random matrices, incidence matrices of undirected graphs, matrices for S-path problem and incidence matrices of directed graphs.

## rationalNumberField/WeightedLinearMatroidParity.cpp 
the algorithm over rational number field

## GaloisField/WeightedLinearMatroidParity.cpp
the algorithm over Galois Field
Also it supports partial pairng (Parts of columns are paired).

## Test cases
Each test case starts with the number of rows and the number of columns, and elements of matrix follow.

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

