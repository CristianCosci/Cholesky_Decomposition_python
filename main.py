import tester as Tester
import argparse
import logging


def main(args):
    # impostazioni per logging
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=logging.INFO if args.verbose else logging.ERROR
        )

    # impostazioni per il test
    Tester.set_algorithm(args.algorithm)

    # questo funziona solo da python 3.10 in poi
    # TODO: ci piace ??
    match args.test_mode:
        case "find_limit":
            Tester.find_limit(seed=args.seed, method=args.method, jit=args.jit)

        case "simple":
            _ = Tester.simple_test(
                    size=args.size, 
                    seed=args.seed, 
                    method=args.method, 
                    jit=args.jit
                )
        case "benchmark":
            _ = Tester.benchmark(
                    size=args.size,
                    seed=args.seed,
                    method=args.method,
                    jit=args.jit
                )

        case _:
            return -1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-tm",
        "--test_mode", 
        type=str,
        choices=["simple", "find_limit", "benchmark"],
        default="simple",
        help="Avvia una tra le varie procedure di testing."
    )

    parser.add_argument(
        "-m",
        "--method", 
        type=str,
        choices=["row", "column", "diagonal"],
        default="column",
        help="Indica l'implementazione di Cholesky da utilizzare."
    )

    parser.add_argument(
        "--jit", 
        action="store_true",
        help="Utilizza la JIT Compile per migliorare le performance."
    )

    parser.add_argument(
        "--seed", 
        type=int,
        default=20,
        help="Specifica il seed per la ripetibilità dei test."
    )

    parser.add_argument(
        "--size", 
        type=int,
        default=10_000,
        help="Specifica la dimensione della matrice (se possibile)."
    )

    parser.add_argument(
        "-alg",
        "--algorithm", 
        type=str,
        choices=["cholesky", "gauss"],
        default="cholesky",
        help="Indica quale tra i possibili algoritmi utilizzare."
    )

    parser.add_argument(
        "-v",
        "--verbose", 
        action="store_true",
        help="Abilita la modalità verbose mostrano molte più informazioni."
    )

    args = parser.parse_args()

    main(args)
