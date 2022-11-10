from time import time
from typing import Any, Callable, Tuple
import numpy as np


def __get_time(millis=True) -> int:
    actual_time = round(time())
    return (actual_time * 1000) if millis else actual_time


def get_execution_time(function: Callable, parameters=[]) -> Tuple[int, Any]:
    start_time = __get_time() # in millisecondi

    result = function(*parameters)

    end_time = __get_time()

    execution_time = end_time - start_time

    return (execution_time, result)


if __name__ == "__main__":
    from time import sleep

    def test(a, b):
        print(f"a: {a}; b: {b}")
        sleep(10)
        return 100

    execution_time, result = get_execution_time(test, [1, 2])

    print(execution_time, result)