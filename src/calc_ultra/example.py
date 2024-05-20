"""
This example is meant to illustrate how to import Calc-Ultra.
In real practice, only line 7 is needed.
"""

if __name__ == "__main__":
    try:
        from calc_ultra import main

        print("\nSuccess.")

    except ImportError:
        print("\nFailed.\nCalc-Ultra is not installed.\n")
