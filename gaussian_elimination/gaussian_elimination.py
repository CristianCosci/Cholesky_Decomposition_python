import numpy as np
from tqdm import tqdm


def compute(Ab: np.ndarray) -> np.ndarray:
        '''
            Gauss elmination algorithm

            inputs:
                Ab ->   Augmented Matrix, matrix of coefficients with the column 
                        of known terms added.

            return:
                The Upper Triangular matrix obtained by applying the Gaussian
                elimination to the given matrix Ab.
            
            N.B.: the returning matrix is not the same as Cholesky, 
                so the operation UU '= A is not valid !!
        '''
        # make a (true) copy of the array
        # and casts it to float
        matrix = Ab.copy()
     
        if matrix[0,0] == 0.0:
            return None
            # raise Exception("matrix row 1 column 1 cannot be zero!")

        n,m = matrix.shape
        
        for i in tqdm(range(0,n), "Gaussian Elimination"):#row
            for j in range(i+1,n):
                if matrix[j,i] != 0.0:
                    matrix[j,i:m]=matrix[j,i:m] - (matrix[j,i]/matrix[i,i])*matrix[i,i:m]

        return matrix


def is_correct_solution(A: np.ndarray, G_U: np.ndarray, b: np.array) -> bool:
    '''
        Check that the result of the Gaussian elimination is correct.

        To do this, the system is solved and verified that 
        the solution obtained is indeed a correct solution to the system.

               ??
            Ax == b

        input:
            A -> Coefficients Matrix
          G_U -> Matrix obtained by applying the Gaussian elimination
            b -> Known System Terms

        return:
            True, if the solution is correct. False otherwise
    '''
    n, _ = G_U.shape
    x = np.zeros(n)

    # calculates the solution of the system with forward substitution
    x[n-1] = G_U[n-1][n]/G_U[n-1][n-1]

    for i in tqdm(range(n-2,-1,-1), "Checking Gauss Solution"):
        x[i] = G_U[i][n]
        
        for j in range(i+1,n):
            x[i] = x[i] - G_U[i][j]*x[j]
        
        x[i] = x[i]/G_U[i][i]

    b_bis = A.dot(x)

    return np.allclose(b_bis, b, 0.001, 0.001)


if __name__ == "__main__":
    # all this server to import (violently)
    # the tester module
    import os
    import sys
    import inspect
    currentdir = os.path.dirname(
        os.path.abspath(
            inspect.getfile(
                inspect.currentframe()
                )
            )
        )
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir) 

    import utils as Tester

    A, b = Tester.generate_data(size=10, seed=20)
    Ab = np.c_[A, b]

    U = compute(Ab)

    print(is_correct_solution(A, U, b))