from typing import Tuple
import numpy as np
import logging


__seed_setted = False


def __generate_A(size=10, seed=None) -> np.ndarray:
    '''
        Genera una matrice Quadrata, Simmetrica e Definita Positiva di dimensione
        size.
    '''
    logging.info(f"Generating A: {size}x{size}")

    # magic ✨
    A = np.random.rand(size, size)
    B = np.dot(A, A.transpose())

    return B


def __generate_b(size=10, seed=None) -> np.array:
    '''
        Genera il vettore dei termini noti
    '''
    
    logging.info(f"Generating b: {size}")

    b = np.random.rand(size)

    return b 


def generate_data(size=10, seed=None) -> Tuple[np.ndarray, np.array]:
    '''
        Ritorna una tupla con i dati del sistema lineare A e b.
    '''
    
    global __seed_setted

    # se esplicitamente passato, viene settato il seed
    if seed is not None and not __seed_setted:
        np.random.seed(seed)

    A =  __generate_A(size, seed)
    b = __generate_b(size, seed)

    __seed_setted = True    # non setta più il seed in futuro
    
    return (A, b)