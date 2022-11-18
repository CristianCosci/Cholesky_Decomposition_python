import utils as Tester
import argparse
import logging

def main(args):
    # settings for logging
    logging.basicConfig(
        format='%(levelname)s: %(message)s', 
        level=logging.INFO if args.verbose else logging.ERROR
        )

    # test settings
    Tester.set_algorithm(args.algorithm)

    # this only works from python 3.10 onwards
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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-tm",
        "--test_mode", 
        type=str,
        choices=["simple", "find_limit", "benchmark"],
        default="simple",
        help=
        """Start the selected test mode.
        simple:     generate data, compute factorization/decomposition and resolve the Linear System.
        find_limit: compute different Cholesky Factorization over bigger matrix (size * 2) every time, starting from a 100x100.
        benchmark:  generate data and only compute the factorization/decomposition. 
                    This returns the execution time and saves results in a file

        """
    )

    parser.add_argument(
        "-m",
        "--method", 
        type=str,
        choices=["row", "column", "diagonal"],
        default="column",
        help="Select which Cholesky implementation to use."
    )

    parser.add_argument(
        "--jit", 
        action="store_true",
        help="Enable JIT to enhance the performance."
    )

    parser.add_argument(
        "--seed", 
        type=int,
        default=20,
        help="Set the seed for the Random Number Generation."
    )

    parser.add_argument(
        "--size", 
        type=int,
        default=10_000,
        help="Set the matrix size (if possible)."
    )

    parser.add_argument(
        "-alg",
        "--algorithm", 
        type=str,
        choices=["cholesky", "gauss"],
        default="cholesky",
        help="Choose the algorithm to use."
    )

    parser.add_argument(
        "-v",
        "--verbose", 
        action="store_true",
        help="Enable verbose mode."
    )

    args = parser.parse_args()

    main(args)
