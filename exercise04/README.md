# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 04

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
To perform the SCF algorithm, we need to produce a unitary transformation matrix that
orthogonalizes the basis functions in the basis set. We will use the canonical procedure
for this which is demonstrated in this exercise.

The unitary transformation matrix can be readily produced by diagonalizing the
overlap matrix and using the inverse of the square root of its eigenvalues as well
as its eigenvectors. Matrix diagonalization is readily available in Python via
the `numpy` library. Note that because we know that our matrices will be Hermitian
(to be more precies, since we only have real-valued elements they will be symmetric),
we can make use of the faster `eigh` as opposed to the more general `eig` function.

## Matrix diagonalization

Performing a matrix diagonalization is then done using
```python
e,v = np.linalg.eigh(M)
```
where `e` is a array of eigenvalues and `v` is a matrix of eigenvectors as columns
by which the following should be equal to the original matrix M

```python
v.transpose().dot(np.diag(e)).dot(v)
```

***Important: Do not use `*` to multiple two matrices. `*` corresponds to an element-wise multiplication which is not what you want. Use `.dot` instead.***

## Question 1

Perform a matrix diagonalization on the overlap matrix for the H2 molecule as constructed
in exercise 3. Show that the above-mentioned equality holds.

## Unitary transformation matrix

The canonical unitary transformation matrix is produced via the following equation

```python
X = v.dot(np.diag(1.0 / np.sqrt(e)))
```

Produce the (canonical) unitary transformation matrix and proof that when this matrix is applied
to the overlap matrix, a unit matrix obtained. 

```python
Sprime = X.transpose().dot(S).dot(X)
```

The way to interpret this result is as follows. The unitary transformation matrix 
orthogonalizes the basis set. This means that a new basis set is produced by means 
of defining a new set of basis functions. This is done by taking for each new basis 
function a linear combination of the known basis functions. In other words, each 
coefficient in the unitary transformation matrix X tells us the linear scaling coefficient 
of the old basis function in the new basis set.

## Solution
The solution is given in [ex04_solution.py](ex04_solution.py)
