from numba import jit
import numpy as np
from tqdm import tqdm
import logging


def compute(A: np.ndarray, method="column", jit=False, nocontrols=False) -> np.ndarray:
    '''
        Apply Cholesky's factoring to obtain the L matrix.
        In order to have a correct computation, the conditions imposed by the
        previous functions have to be respected.

        jit: if True, optimize the execution of this function with numba.
    '''
    # if specified, the initial checks are skipped to save time
    if nocontrols: 
        logging.info("Skipping Requirements")
        is_factorizable = True
    
    else:
        is_factorizable = __check_requirements(A)

    # the constraints are not satisfied, the given matrix cannot be factored
    if not is_factorizable:
        return None

    logging.info(f"Computing Cholesky Factorization {method} - jit: {jit}")
    L = methods[method](A, jit) # avvia la relativa implementazione

    return L


def is_correct_solution(A: np.ndarray, L: np.ndarray) -> bool:
		'''
            Check that the solution is correct by recalculating A from L
		'''
		A_bis = np.dot(L, np.transpose(L))
        # print(A_bis)

		return np.allclose(A, A_bis, 0.001, 0.001)


def __check_requirements(A: np.ndarray) -> bool:

    def is_square(matrix: np.ndarray) -> bool:
        '''
            This function checks that the given matrix is Square:
                A = nxn
        '''
        logging.info("Checking IS SQUARE")
        n, m = matrix.shape
        if(m == n):
            return True
        
        return False

    def is_symmetric(matrix: np.ndarray) -> bool:
        '''
            This function checks that the matrix is Symmetric:

                a_ij = a_ji     (must be symmetrical on the diagonal)

                    oppure
                
                At = A          (the transpose is equal to itself)
        '''
        logging.info("Checking IS SYMMETRIC")

        if(np.allclose(matrix.transpose(), matrix)):
            return True
        
        return False

    def is_positive_definite(matrix: np.ndarray) -> bool:
        '''
            This function checks that the matrix is Definite Positive:

                eigenvalues > 0     (all eigenvalues of the matrix must be positive)
        '''
        logging.info("Checking IS POSITIVE DEFINITE")

        eigenvals = np.linalg.eigvals(matrix)

        if (np.all(eigenvals > 0)):
            return True

        return False


    logging.info("Checking Cholesky Requirements")
    
    # check if all the requirements are met
    return is_square(A) and is_symmetric(A) and is_positive_definite(A)


# --- JIT COMPILED FUNCTIONS --- #

# Note: it seems that the functions that are subject to the JIT cannot fit 
# inside the function that does the factoring since this makes the execution 
# performance worse (on my pc by 1000 ms).
#
# TODO: controllare che effettivamente è vero e il perchè.
# TODO: nelle seguenti funzioni compilate abbiamo una roba del tipo
#       L[i, j] = ...
#       Sta cosa funziona perchè L è passato per riferimento, ma ci piace come cosa ?

## ~~ by COLUMN
@jit(nopython=True)
def __column_columns_numba(L: np.array,A: np.array,i: int,j: int) -> bool:
    L[i,j] = np.sqrt(A[i,j]-np.sum(L[i,:j]**2))

@jit(nopython=True)
def __row_columns_numba(L: np.array,A: np.array,i: int,j: int) -> bool:
    L[i,j] = (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]


## ~~ by ROW
# Calculate the factorization of the elements along the diagonal of A
@jit(nopython=True)
def __row_diagonal(L:np.ndarray, A:np.ndarray, i:int, j:int) -> bool:
    L[i,j] = np.sqrt(A[i,j]-np.sum(L[:i,j]**2))

# Calculate the factorization of elements that are not on the diagonal of A
@jit(nopython=True)
def __row_non_diagonal(L:np.ndarray, A:np.ndarray, i:int, j:int) -> bool:
    L[i,j] = (A[i,j]-np.sum(L[:i,j]*L[:i,i])) / L[i,i]


## ~~ by DIAGONAL
@jit(nopython=True)
def __cholesky_formula_diagonal(i, j, A, L):
        if (i == j):
            return np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calculate the values of the diagonal
        
        return (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]    # calculate the column values


# --- CHOLESKY METHODS --- #

def __compute_by_column(A: np.ndarray, jit=False) -> np.ndarray:
    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)    # initialize the result matrix

    # this if is ugly but maybe it is necessary. 
    # It could be put inside the for but in this way 
    # i would go to perform a check many times that must be performed only once (at the beginning).
    if jit:
        for j in tqdm(range(n), "Cholesky - COLUMN (JIT)"):
            for i in range(j, n):
                if (i == j):
                    __column_columns_numba(L,A,i,j)   # calculate the values of the diagonal
                else:
                    __row_columns_numba(L,A,i,j)   # calculate the column values
    
    else:
        for j in tqdm(range(n), "Cholesky - COLUMN"):
            for i in range(j, n):
                if (i == j):
                    L[i,j] = np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calculate the values of the diagonal
                else:
                    L[i,j] = (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]   # calculate the column values

    return L


def __compute_by_row(A: np.ndarray, jit=False) -> np.ndarray:
    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)  # initialize the result matrix

    if jit:
        for i in tqdm(range(n), "Cholesky - ROW (JIT)"):
            for j in range(i, n):
                if (i == j):
                    __row_diagonal(L,A,i,j)   # calculate the values of the diagonal
                else:
                    __row_non_diagonal(L,A,i,j)  # calculate the column values

    else:
        for i in tqdm(range(n), "Cholesky - ROW"):
            for j in range(i, n):
                if (i == j):
                    L[i,j] = np.sqrt(A[i,j]-np.sum(L[:i,j]**2))   # calculate the values of the diagonal
                else:
                    L[i,j] = (A[i,j]-np.sum(L[:i,j]*L[:i,i])) / L[i,i]   # calculate the column values

    # N.B. the matrix must be transposed !!
    return L.transpose()


def __compute_by_diagonal(A: np.ndarray, jit=False) -> np.ndarray:
    # TODO: is what works but could be improved (From an aesthetic point of view)!

    def cholesky_formula(i, j, A, L):
        if (i == j):
            return np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calculate the values of the diagonal
        
        return (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]    # calculate the column values


    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)    # initialize the result matrix

    external = 0 #variable to count how many external loops I have to do
    internal = 2 * n - 1 - 1
    aux = 0
    if jit:
        for row in tqdm(range(2 * n - 1), "Cholesky - DIAGONAL (JIT)"):
            if row < n-1:
                col = 0
                L[row, col] = __cholesky_formula_diagonal(row, col, A, L)
                aux += 1
                for z in range(1, int(np.floor(row/2))+1):
                    L[row-z, col+z] = __cholesky_formula_diagonal(row-z, col+z, A, L)

            else:
                col = n-1
                L[col, row-aux] = __cholesky_formula_diagonal(col, row-aux, A, L)
                for z in range(1, int(np.floor(internal/2))+1):
                    L[col-z, row-aux+z] = __cholesky_formula_diagonal(col-z, row-aux+z, A, L)

            internal -= 1   
            external += 1
    
    else:
        for row in tqdm(range(2 * n - 1), "Cholesky - DIAGONAL"):
            if row < n-1:
                col = 0
                L[row, col] = cholesky_formula(row, col, A, L)
                aux += 1
                for z in range(1, int(np.floor(row/2))+1):
                    L[row-z, col+z] = cholesky_formula(row-z, col+z, A, L)

            else:
                col = n-1
                L[col, row-aux] = cholesky_formula(col, row-aux, A, L)
                for z in range(1, int(np.floor(internal/2))+1):
                    L[col-z, row-aux+z] = cholesky_formula(col-z, row-aux+z, A, L)

            internal -= 1   
            external += 1
        
    return L


methods = {"row": __compute_by_row, "column": __compute_by_column, "diagonal": __compute_by_diagonal}
