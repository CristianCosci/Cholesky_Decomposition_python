from asyncio.log import logger
import json
from typing import Any, Callable, Dict, Tuple
from .execution_time import get_execution_time
from .data_generator import generate_data
import cholesky_factorization as Cholesky_factorization
import gaussian_elimination as Gaussian_elimination
import linear_system_solver as Linear_sistem
import numpy as np


ALGORITHM = "cholesky"


def simple_test(size=100, seed=20, method="column", jit=False) -> Tuple[np.array, Tuple[int, int]]:
    '''
        Risolve un dato Sistema Lineare con la fattorizzazione di Cholesky.

        inputs:
            A:          matrice che rappresenta il Sistema
            b:          vettore dei termini noti
            jit:        applica la JIT per migliorare le performance di Cholesky

        returns:
            Tuple(
                x -> soluzione del sistema
                Tuple (
                    c,  -> tempo (in ms) impiegato per risolvere Cholesky
                    l   -> tempo (in ms) impiegato per risolvere il Sistema
                )
            )
    '''
    
    # print(f"fun: {ALGORITHM}")

    def cholesky():
        # --- Visualizzo i dati iniziali --- #
        print(f"A:\n{A}\n")
        print(f"b:\n{b}\n")

        # --- Decomposizione della matrice A in LU --- #
        print("CHOLESKY FACTORIZATION ...")

        cholesky_execution_time, L = get_execution_time(Cholesky_factorization.compute, [A, method, jit])
        
        if L is None:
            print("Impossibile scomporre la matrice data !!", force=True)
            return -1
        
        print(f"L:\n{L}\n")
        print(f"Il risultato è corretto ?: {'✅' if (Cholesky_factorization.is_correct_solution(A, L)) else '❌'}")
        print(f"Tempo di esecuzione: {cholesky_execution_time}")
        print('\n')

        # --- Risoluzione del Sistema Lineare --- #
        print("SOLVING LINEAR SYSTEM ...")

        linsys_execution_time, x = get_execution_time(Linear_sistem.solve, [L, None, b])

        print(f"x: \n{x}\n")
        print(f"Il risultato è corretto ?: {'✅' if (Linear_sistem.is_correct_solution(A, x, b)) else '❌'}")
        print(f"Tempo di esecuzione: {linsys_execution_time}")

        return (x, (cholesky_execution_time, linsys_execution_time))

    def gauss():
        Ab = np.c_[A, b]    # Augmented Matrix

        # --- Visualizzo i dati iniziali --- #
        print(f"A:\n{A}\n")
        print(f"b:\n{b}\n")
        print(f"Ab:\n{Ab}\n")

        # --- Applicazione dell'algoritmo di eliminaizone di Gauss --- #
        print("GAUSSIAN ELIMINATION ...")

        gauss_execution_time, U = get_execution_time(Gaussian_elimination.compute, [Ab])
        
        if U is None:
            print("Impossibile applicare l'algoritmo di Gauss sualla matrice data !!")
            return -1
        
        print(f"U:\n{U}\n")
        print(f"Il risultato è corretto ?: {'✅' if (Gaussian_elimination.is_correct_solution(A, U, b)) else '❌'}")
        print(f"Tempo di esecuzione: {gauss_execution_time}")
        print('\n')

        # --- Risoluzione del Sistema Lineare --- #
        print("SOLVING LINEAR SYSTEM ...")

        linsys_execution_time, x = get_execution_time(Linear_sistem.solve, [None, U, None])

        print(f"x: \n{x}\n")
        print(f"Il risultato è corretto ?: {'✅' if (Linear_sistem.is_correct_solution(A, x, b)) else '❌'}")
        print(f"Tempo di esecuzione: {linsys_execution_time}")

        return (x, (gauss_execution_time, linsys_execution_time))


    A, b = generate_data(size=size, seed=seed)

    match ALGORITHM:
        case "cholesky":
            return cholesky()
        case "gauss":
            return gauss()
        case _:
            raise Exception("Bad Algorithm Name !")


def find_limit(starting_size=100, seed=20, method="column", jit=False):
    size = starting_size
    
    while True:
        print(f"SIZE: {size}x{size}")

        A, b = generate_data(size, seed)

        execution_time, _ =  get_execution_time(Cholesky_factorization.compute, [A, method, jit])
        
        print(f"Cholesky Execution Time: {execution_time} ms")
        print("\n")
        
        size *= 2


def benchmark(size=10_000, seed=20, method="column", jit=False) -> Tuple[int, Any]:  # TODO: controllare cosa ritornare
    print("Generating data ...")
    A, b = generate_data(size, seed)
    print()

    match ALGORITHM:
        case "cholesky":
            print(f"ALGORITHM:\t Cholesky (by {method})")
            print(f"MATRIX SIZE:\t {size}x{size}")
            print(f"SEED:\t\t {seed}")
            print(f"JIT:\t\t {jit}")
            print("")

            execution_time, L =  get_execution_time(Cholesky_factorization.compute, [A, method, jit, True])
            print(f"Execution Time: {execution_time}")

            data = {
                #'A': A.tolist(), 
                #'b': b.tolist(),
                #'res': L.tolist(),
                'time': execution_time,
                'algorithm': ALGORITHM, 
                'size': size, 'seed': seed, 
                'method': method,
                'jit': jit
                }

            __save(data)

        case "gauss":
            print(f"ALGORITHM:\t Gauss")
            print(f"MATRIX SIZE:\t {size}")
            print(f"SEED:\t\t {seed}")
            print("")

            Ab = np.c_[A, b]    # Augmented Matrix
            execution_time, U = get_execution_time(Gaussian_elimination.compute, [Ab])
            print(f"Execution Time: {execution_time}")

            data = {
                #'A': A.tolist(), 
                #'b': b.tolist(),
                #'res': U.tolist(),
                'time': execution_time,
                'algorithm': ALGORITHM, 
                'size': size, 'seed': seed, 
                'method': None,
                'jit': None,
                }

            __save(data)

        case _:
            raise Exception("Bad Algorithm Name !")


def set_algorithm(string: str):
    '''
        Cambia l'algoritmo da utilizzare

            Cholesky/Gauss
    '''
    global ALGORITHM

    ALGORITHM = string


def __save(data: Dict):
    logger.info("Saving Data")
    with open(f"{data['algorithm']}_{data['method']}_{data['size']}_{data['seed']}.json", "w") as f:
        jsonobj = json.dumps(data)
        f.write(jsonobj)