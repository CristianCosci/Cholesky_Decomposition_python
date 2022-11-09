import numpy as np
from tqdm import tqdm


def compute(Ab: np.ndarray) -> np.ndarray:
        '''
            Gauss elmination algorithm

            inputs:
                Ab ->   Augmented Matrix, matrice dei coefficienti con aggiunta
                        la colonna dei termini noti.

            return:
                La matrice Triangolare Superiore ottenuta applicando l'eliminazione
                Gaussiana alla matrice data Ab.
            
            N.B.: la matrice che ritorna non è la stessa cosa di Cholesky, quindi
                    l'operazione UU' = A non vale !!
        '''

        # effettua una copia (vera) dell'array
        # e lo casta a float
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
        Controlla che il risultato dell'eliminazione gaussiana sia corretto.

        Per farlo viene risolto il sistema e verificato che la soluzione ottenuta
        sia effettivamente una soluzione corretta al sistema.

               ??
            Ax == b

        input:
            A -> Matrice dei Coefficienti
          G_U -> Matrice ottenuta dall'applicazione dell'eleiminazione gaussiana
            b -> Termini noti del Sistema

        return:
            True, se la soluzione è corretta. False altrimenti
    '''

    n, _ = G_U.shape
    x = np.zeros(n)

    # calcola la soluzione del sistema con forward sostitution
    x[n-1] = G_U[n-1][n]/G_U[n-1][n-1]

    for i in tqdm(range(n-2,-1,-1), "Checking Gauss Solution"):
        x[i] = G_U[i][n]
        
        for j in range(i+1,n):
            x[i] = x[i] - G_U[i][j]*x[j]
        
        x[i] = x[i]/G_U[i][i]

    b_bis = A.dot(x)

    return np.allclose(b_bis, b, 0.001, 0.001)


if __name__ == "__main__":
    # tutto questo server per importare (in modo violento)
    # il modulo tester
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

    import tester as Tester

    A, b = Tester.generate_data(size=10, seed=20)
    Ab = np.c_[A, b]

    U = compute(Ab)

    print(is_correct_solution(A, U, b))