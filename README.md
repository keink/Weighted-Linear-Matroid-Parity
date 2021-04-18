# Description
This is an implementation of an algorithm for weighted linear matroid parity problem, proposed by Iwata and Kobayashi [1].

# Requirements
[Boost Multiprecision Library](https://www.boost.org/doc/libs/1_75_0/libs/multiprecision/doc/html/index.html)

# Usage
+ You can switch types of matrices by standard input (it is inputted to a variable "inc").
Type 0, 1, 2, 4, 10 respectively correspond to integer random matrices, incidence matrices of undirected graphs, matrices for S-path problem and incidence matrices of directed graphs, and for approximation algorithms for steiner tree problem (only available for the GaloisField program).

+ You can run programs by replacing filepaths in main function with appropriate paths where files of test cases are located.

+ Three variales "elements", "row_size", "nonzero_probs" are used for controlling to what extent the algorithm is applied among all the sample test cases (available in "sotsuron_data" directory).  
Those variables respectively represent:  
  + "elements": absolute value of elements in matrices  
  + "row_size": the number of rows of matrices  
  + "nonzero_prob": the probabilities of nonzero elements (used to create ramdom sample test cases)  
If you want to use your own test cases not sample test cases, then you can removes for-loops on these variables in main function.


## rationalNumberField/WeightedLinearMatroidParity.cpp 
the algorithm over rational number field

## GaloisField/WeightedLinearMatroidParity.cpp
the algorithm over Galois Field
It supports partial pairng (Not all columns are necessarily paired to form lines).

## Test cases
Each test case starts with the number of rows and the number of columns, follwed by elements of matrix and weights of lines.
For instance, if you want to run the algorithm for 
<img src = https://latex.codecogs.com/gif.latex?\begin{pmatrix}&space;0&space;&&space;1&space;&&space;2&space;&&space;3&space;\\&space;-1&space;&&space;0&space;&&space;1&space;&&space;2&space;\\&space;-2&space;&&space;-1&space;&&space;0&space;&&space;1&space;\\&space;-3&space;&&space;-2&space;&&space;-1&space;&&space;0&space;\end{pmatrix} />
where weights of lines are (1, 2), then test case should be written as below.
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

