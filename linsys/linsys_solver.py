import numpy as np
from tqdm import tqdm


def solve(L: np.ndarray=None, G_U: np.ndarray=None, b: np.array=None) -> np.array:
    '''
        Risolve il Sistema Lineare dato.
    
        In base ai parametri che ottine:

            - se viene dato solo G_U, risolve il sistema partendo dalla matrice
                ottenuta con l'eliminazione Gaussiana
            - se viene dato solo L e b, risolve il sistema partendo dalla matrice
                L ottenuta con la fattorizzazione di Cholesky

            - se i parametri non rispettano i critri precedenti,
                lancia un'eccezione
    '''

    if L is not None and b is not None and G_U is None:
        return __solve_cholesky(L, b)

    if G_U is not None and L is None and b is None:
        return __solve_gauss(G_U)

    raise Exception("Wrong Parameters")


def is_correct_solution(A: np.ndarray, x: np.array, b: np.array) -> bool:
    '''
        Controlla che la soluzione x al sistema dato Ab sia corretta.
        Questo viene fatto moltiplicando A per x e controllando che il risultato sia
        uguale a b
               
               ??
            Ax == b
    '''

    b_bis = A.dot(x)

    return np.allclose(b, b_bis, 0.001, 0.001)    # TODO: testare


def __solve_cholesky(L: np.ndarray, b: np.array) -> np.array:
    '''
        Solve linear sistem using cholesky decomposition 
        with matrix L
        
        L -> Decomposition matrix by cholesky ( lower triangular matrix ) 
        U -> Upper triangular matrix 
        b -> Know term 
        x,y -> The solutio of equation
        x -> The solution using Backword subsostitution
        y -> The solution using Forword sostitution
        
        a_11 a12 ... a_1n  x_1      b1
        a_12 a22 ... a_2n  x_2      b2
        ...  ... ... ...         =  ...
        a_n1 an2 ... a_nn  x_n      bn 
        
        linear system --->   [a] {x} = {b}
        
        1. Decomposition: [a] = [L][U]
        2. Substitution: [L][U] = {b}
            2.1 Backword substitution [U]{x} = {y}
            2.2 Forward substitution [L]{y} = {b}  

        THE SYSTEM IS SOLVABLE? 

            In the case where A is full rank, the matrix A^TA is symmetrical defined positive, 
            and it is therefore possible to solve the normal system by Cholesky factorization.

            In the case where A is not symmetrical defined positive is not decomposable by Cholesky
            and the linear system can be solved in other ways but not with Cholesky. 
               
    '''

    U = np.transpose(L)
    n, _ = np.shape(L)

    # soluzioni (forword/backword)
    y = np.zeros(n, dtype=np.float64)   # forzo tipo float64
    x = np.zeros(n, dtype=np.float64)

    # TODO: si potrebbe aggiungere un argomento per scegliere quale 
    #       metodo utilizzare.

    # Forword sostitution
    for i in tqdm(range(n), "Solving Cholesky (Forword)"):
        sumj = 0                    # TODO: si portrebbe riscrivere con la 
        for j in range(i):          # sommatoria di numpy
            sumj += L[i, j] * y[j]

        y[i] = (b[i]-sumj)/L[i, i]

    # Backword sostitution
    for i in tqdm(range(n-1, -1, -1), "Solving Cholesky (Backword)"):
        sumj = 0                    # TODO: si portrebbe riscrivere con la 
        for j in range(i+1, n):     # sommatoria di numpy
            sumj += U[i, j] * x[j]

        x[i] = (y[i]-sumj)/U[i, i]

    return x
    # return (x, y)


def __solve_gauss(G_U: np.ndarray) -> np.array:
    n, _ = G_U.shape
    x = np.zeros(n)

    # calcola la soluzione del sistema con forward sostitution
    x[n-1] = G_U[n-1][n]/G_U[n-1][n-1]

    for i in tqdm(range(n-2,-1,-1), "Solving Gauss"):
        x[i] = G_U[i][n]
        
        for j in range(i+1,n):
            x[i] = x[i] - G_U[i][j]*x[j]
        
        x[i] = x[i]/G_U[i][i]

    return x