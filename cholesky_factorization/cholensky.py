from numba import jit
import numpy as np
from tqdm import tqdm
import logging


def compute(A: np.ndarray, method="column", jit=False, nocontrols=False) -> np.ndarray:
    '''
        Applica la fattorizzazinoe di Cholesky per ottenere la matrice L
        Per poter funzionare devono essere rispettate le condizioni imposte dalle
        precedenti funzioni.

        jit:   se True, ottimizza l'esecuzione di questa funzione con numba.
    '''

    # se specificato i controlli iniziali vengono saltati
    # per risparmiare tempo
    if nocontrols: 
        logging.info("Skipping Requirements")
        is_factorizable = True
    
    else:
        is_factorizable = __check_requirements(A)

    # i vincoli non sono soddisfatti, la matrice data non si può fattorizzare
    if not is_factorizable:
        return None

    logging.info(f"Computing Cholesky Factorization {method} - jit: {jit}")
    L = methods[method](A, jit) # avvia la relativa implementazione

    return L


def is_correct_solution(A: np.ndarray, L: np.ndarray) -> bool:
		'''
            Controlla che la soluzione sia corretta ricalcolando A da L
		'''

		A_bis = np.dot(L, np.transpose(L))
		
        # print(A_bis)

		return np.allclose(A, A_bis, 0.001, 0.001)


def __check_requirements(A: np.ndarray) -> bool:

    def is_square(matrix: np.ndarray) -> bool:
        '''
            Questa funzione controlla che la matrice data sa Quadrata:

                A = nxn
        '''

        logging.info("Checking IS SQUARE")

        n, m = matrix.shape

        if(m == n):
            return True
        
        return False

    def is_symmetric(matrix: np.ndarray) -> bool:
        '''
            Questa funzione controlla che la matrice si Simmetrica:

                a_ij = a_ji     (deve essere simmetrica sulla diagonale)

                    oppure
                
                At = A          (la trasposta è uguale a se stessa)

                
        '''

        logging.info("Checking IS SYMMETRIC")


        if(np.allclose(matrix.transpose(), matrix)):
            return True
        
        return False

    def is_positive_definite(matrix: np.ndarray) -> bool:
        '''
            Questa funzione controlla che la matrice sia Definita Positiva:

                eigenvalues > 0     (tutti gli eigenvalue della matrice devono essere positivi)
        '''

        logging.info("Checking IS POSITIVE DEFINITE")


        eigenvals = np.linalg.eigvals(matrix)

        if (np.all(eigenvals > 0)):
            return True

        return False


    logging.info("Checking Cholesky Requirements")
    
    # controlla se tutti i requisiti sono soddisfatti
    return is_square(A) and is_symmetric(A) and is_positive_definite(A)


# --- JIT COMPILED FUNCTIONS --- #

# Note: sembra che le funzioni che sono soggette alla JIT non possono stare
#       all'interno della funzione che effettua la fattorizzazione dato che
#       ciò fa peggiorare le performance di esecuzione (sul mio pc di 1000 ms).
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
# Calcola la fattorizzazione degli elementi lungo la diagonale di A
@jit(nopython=True)
def __row_diagonal(L:np.ndarray, A:np.ndarray, i:int, j:int) -> bool:
    L[i,j] = np.sqrt(A[i,j]-np.sum(L[:i,j]**2))

# Calcola la fattorizzazione degli elementi che non sono sulla diagonale di A
@jit(nopython=True)
def __row_non_diagonal(L:np.ndarray, A:np.ndarray, i:int, j:int) -> bool:
    L[i,j] = (A[i,j]-np.sum(L[:i,j]*L[:i,i])) / L[i,i]


## ~~ by DIAGONAL
@jit(nopython=True)
def __cholesky_formula_diagonal(i, j, A, L):
        if (i == j):
            return np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calcolo i valori della diagonale
        
        return (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]    # calcolo i valori delle colonne


# --- CHOLESKY METHODS --- #

def __compute_by_column(A: np.ndarray, jit=False) -> np.ndarray:
    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)    # inizializzo la matricce risultato

    # questo if è bruttino ma forse è necessario. Potrebbe essere messo all'interno
    # del for ma così andrei ad effettuare tante volte un controllo che deve 
    # essere eseguito una volta solo (all'inizio).
    if jit:
        for j in tqdm(range(n), "Cholesky - COLUMN (JIT)"):
            for i in range(j, n):
                if (i == j):
                    __column_columns_numba(L,A,i,j)   # calcolo i valori della diagonale
                else:
                    __row_columns_numba(L,A,i,j)   # calcolo i valori delle colonne
    
    else:
        for j in tqdm(range(n), "Cholesky - COLUMN"):
            for i in range(j, n):
                if (i == j):
                    L[i,j] = np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calcolo i valori della diagonale
                else:
                    L[i,j] = (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]   # calcolo i valori delle colonne

    return L


def __compute_by_row(A: np.ndarray, jit=False) -> np.ndarray:
    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)  # inizializzo la matricce risultato

    if jit:
        for i in tqdm(range(n), "Cholesky - ROW (JIT)"):
            for j in range(i, n):
                if (i == j):
                    __row_diagonal(L,A,i,j)   # calcolo i valori della diagonale
                else:
                    __row_non_diagonal(L,A,i,j)  # calcolo i valori delle colonne

    else:
        for i in tqdm(range(n), "Cholesky - ROW"):
            for j in range(i, n):
                if (i == j):
                    L[i,j] = np.sqrt(A[i,j]-np.sum(L[:i,j]**2))   # calcolo i valori della diagonale
                else:
                    L[i,j] = (A[i,j]-np.sum(L[:i,j]*L[:i,i])) / L[i,i]   # calcolo i valori delle colonne

    # N.B. la matrice va trasposta !!
    return L.transpose()


def __compute_by_diagonal(A: np.ndarray, jit=False) -> np.ndarray:
    # TODO: sta cosa funziona ma potrebbe essere migliorata !

    def cholesky_formula(i, j, A, L):
        if (i == j):
            return np.sqrt(A[i,j]-np.sum(L[i,:j]**2))   # calcolo i valori della diagonale
        
        return (A[i,j]-np.sum(L[i,:j]*L[j,:j])) / L[j,j]    # calcolo i valori delle colonne


    n, _ = A.shape

    L = np.zeros(n*n, dtype=float).reshape(n, n)    # inizializzo la matricce risultato

    external = 0 #variabile per contare quanti cicli esterni devo fare
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
