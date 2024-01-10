exclude = ["pi", "e"]

from sympy.core.numbers import pi, E, oo
from math import floor, ceil, erf, fabs, factorial
from numpy import *

for name in exclude:
    del globals()[name]

from sympy import (
    diff,
    idiff,
    Integral,
    integrate,
    limit,
    pprint,
    simplify,
    symbols,
    pi,
    E,
)
import matplotlib.pyplot as plt
import datetime, logging, os, time, warnings


# This program requires the Sympy, NumPy, and MatPlotLib modules installed!

# Disable auto Python warnings

warnings.filterwarnings("ignore")


def derive(function, order):
    function = replace_expr(function)

    if check_order(order) is True:
        global df
        df = diff(function, x, order)

        nprint(f"\nDerivative of {function} with order {order} is:\n")
        print_expr(df)

        if check_simp(df) is True:
            nprint("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))

        df = trig_rep(str(df))

        nprint("\nShow graph of area? (y/n)")
        show = input("(Exit the graph window when you are finished to continue) ")

        if show == "y":
            try:
                nprint("\nLoading graph. Might take some time on first startup ...")
                x_array = linspace(-8, 8, 200000)

                title = "Function (red) and derivative (blue)"
                plt.title(title)
                plt.xlabel("x", weight="bold")
                plt.ylabel("y", rotation=0, weight="bold")
                plt.plot(x_array, f(x_array), color="red", label="Function")
                plt.plot(x_array, dif(x_array), color="blue", label="Derivative")
                plt.axis([-7.5, 7.5, -7.5, 7.5])
                plt.legend(loc="lower left")
                plt.grid()
                plt.show()

                return "\nExited graph."

            except:
                plt.close()
                logging.warning("Could not graph function.")
                return "\nExited graph."


def partial_derive(function, var, order):
    function = replace_expr(function)

    if check_order(order) is True:
        df = diff(function, var, order)

        nprint(
            f"\nPartial derivative of {function} in respect to {var} of order {order} is:\n"
        )
        print_expr(df)

        if check_simp(df) is True:
            nprint("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))


def implicit_derive(circ, order):
    circ = replace_expr(circ)

    if check_order(order) is True:
        df = idiff(eval(circ), y, x, int(order))

        nprint(f"\nDerivative of {circ} with order {order} is:\n")
        print_expr(df)

        if str(simplify(df, evaluate=False)) != str(df):
            nprint("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))


def antiderive(function):
    function = replace_expr(function)

    global F
    F = Integral(function, x).doit()

    if "Integral" in str(F):
        logging.warning("Cannot compute integral.\n")
        return ""

    nprint("\nAntiderivative is:\n")
    print_expr(F)

    if check_simp(F) is True:
        nprint("\nSimplify/rewrite:\n")
        print_expr(simplify(F, evaluate=False))

    F = trig_rep(str(F))

    nprint("\nDon't forget to add a constant!\n")

    nprint("\nShow graph of area? (y/n)")
    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            nprint("\nLoading graph. Might take some time on first startup ...")

            x_array = linspace(-8, 8, 200000)

            title = "Function (red) and antiderivative (blue, C = 0)"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_array, f(x_array), color="red", label="Function")
            plt.plot(x_array, af(x_array), color="blue", label="Antiderivative")
            plt.axis([-7.5, 7.5, -7.5, 7.5])
            plt.legend(loc="lower left")

            plt.grid()
            plt.show()

            return "\nExited graph."

        except:
            plt.close()
            logging.warning("Could not graph function.")
            return "\nExited graph."


def definite_integrate(function, low, up):
    x = symbols("x")
    function = replace_expr(function)

    if check_bound(low) is False:
        return ""

    low = eval(low)

    if check_bound(up) is False:
        return ""

    up = eval(up)

    result = integrate(function, (x, low, up))

    up = eval(replace_bound(str(up)))
    low = eval(replace_bound(str(low)))

    num_result = integrate(function, (x, low, up)).evalf()

    if (
        (str(result) == "nan")
        or ("I" in str(result))
        and ("Integral" not in str(result))
    ):
        logging.warning("Cannot compute integral because integral does not converge.")
        return ""

    if "Integral" not in str(result):
        nprint(
            f"\nCalculated integral of {function} from {low} to {up}. Final area is:\n"
        )
        print_expr(result)

    else:
        nprint("\nCannot express result symbolically.")

    nprint("\nNumeral evaluation/approximation:\n")
    print_expr(num_result)

    nprint("\nShow graph of area? (y/n)")

    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            nprint("\nLoading graph. Might take some time on first startup ...")

            x_array = linspace((-up - 8), (up + 8), 200000)

            if "log" in function:
                x_array = linspace(
                    int(floor(low)) - 3,
                    int(ceil(up)) + 8,
                    200000,
                )

            title = "Shaded area beneath function"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_array, f(x_array), color="red", label="Function")
            plt.fill_between(
                x_array,
                f(x_array),
                where=[(x_array > low) and (x_array < up) for x_array in x_array],
                color="blue",
            )

            try:
                if graph_option == "f":
                    plt.axis([-7.5, 7.5, -7.5, 7.5])

                elif graph_option == "a":
                    if (float(f(low)) != 0) and (float(f(up)) != 0):
                        plt.axis(
                            [
                                low - 5,
                                up + 5,
                                float(f(round(low)))
                                - (float(f(round(low))) + float(f(round(up)))) / 2
                                - 1,
                                float(f(round(up)))
                                + (float(f(round(low))) + float(f(round(up)))) / 2
                                + 1,
                            ]
                        )

                    elif (float(f(low)) == 0) or (float(f(up)) == 0):
                        plt.axis([low - 5, up + 5, -(up - low) / 2, (up + low) / 2])

            except:
                plt.axis([-7.5, 7.5, -7.5, 7.5])

            plt.legend(loc="lower left")
            plt.grid()
            plt.show()

            return "\nExited graph."

        except:
            plt.close()
            logging.warning("Could not graph function.")
            return "\nExited graph."

    else:
        return "\nExiting Definite Integral Screen ... ... ...\n"


def improper_integrate(function, low, up):
    function = replace_expr(function)

    if check_bound(low) is False:
        return ""

    low = eval(low)

    if check_bound(up) is False:
        return ""

    up = eval(up)

    try:
        improper_area = Integral(function, (x, low, up)).principal_value()

        if "Integral" not in str(improper_area):
            nprint(
                f"\nCalculated improper integral of {function} from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)

        else:
            nprint("Cannot compute improper integral.")

        print()

    except ValueError:
        logging.warning("ValueError: Singularity while computing improper integral.\n")


def normal_limit(expr, value):
    if "pi" in expr:
        expr = expr.replace("pi", str(pi))

    replace_expr(expr)

    if check_bound(value) is False:
        return ""

    value = float(eval(replace_bound(value)))

    l = limit(expr, x, value)

    if "Limit" in str(l):
        logging.warning("Cannot compute limit.")
        return ""

    nprint(f"\nLimit of {expr} as x approaches {value} is:\n")
    print_expr(l)
    if check_simp(l) is True:
        nprint("\nSimplify/rewrite:\n")
        print_expr(simplify(l, evaluate=False))


def one_side_limit(expr, value, direction):
    if "pi" in expr:
        expr = expr.replace("pi", str(pi))

    replace_expr(expr)

    if check_bound(value) is False:
        return ""

    value = float(eval(replace_bound(value)))

    if direction == "left":
        direction_sign = "-"

    elif direction == "right":
        direction_sign = "+"

    else:
        logging.error("\nTypeError: Direction is neither right or left.")
        return ""

    l = limit(expr, x, value, dir=direction_sign)

    if "Limit" in str(l):
        logging.warning("\nCannot compute limit.")
        return ""

    nprint(f"\nLimit of {expr} as x approaches {value} from the {direction} is:\n")
    print_expr(l)


def check_simp(expr):
    if str(simplify(expr, evaluate=False)) != str(expr):
        return True

    else:
        return False


def check_order(order):
    if ("." in order) or (order.isnumeric() == False) or (int(order) <= 0):
        logging.error(
            "OrderError: Order of derivative calculation is not a valid number."
        )
        return False

    else:
        return True


def check_bound(bound):
    if (
        (bound.isnumeric() is False)
        and ("pi" not in bound)
        and ("E" not in bound)
        and ("-" not in bound)
        and ("." not in bound)
        and ("sqrt" not in bound)
        and ("oo" not in bound)
        and ("/" not in bound)
    ):
        logging.error("TypeError: Integration bound is a not a number.")
        return False

    else:
        return True


def replace_expr(expr):
    if "^" in expr:
        expr = expr.replace("^", "**")

    return expr


def replace_bound(bound):
    if "pi" in bound:
        bound = bound.replace("pi", str(pi.evalf()))
    if "E" in bound:
        bound = bound.replace("E", str(E.evalf()))

    return bound


def f(x):
    if "Abs" in dfunction:
        final = dfunction.replace("Abs", "fabs")
    else:
        final = dfunction

    if "x" not in final:
        final = "0 * x + " + final

    final = eval(trig_rep(final))

    return final


def dif(x):
    if "x" not in str(df):
        final = eval("0 * x + " + trig_rep(str(df)))

    else:
        final = eval(df)

    return final


def af(x):
    if "x" not in str(F):
        final = eval("0 * x + " + trig_rep(str(F)))

    else:
        final = eval(F)

    return final


def trig_rep(function):
    if "asin" in function:
        function = function.replace("asin", "arcsin")

    if "acos" in function:
        function = function.replace("acos", "arccos")

    if "atan" in function:
        function = function.replace("atan", "arctan")

    if "asinh" in function:
        function = function.replace("asinh", "arcsinh")

    if "acosh" in function:
        function = function.replace("acosh", "arccosh")

    if "atanh" in function:
        function = function.replace("atanh", "arctanh")

    return function


def print_expr(text):
    printing_methods = {"p": lambda t: pprint(text), "n": lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods["p"](text)


def nprint(text):
    print(text)
    time.sleep(0.04)


def main():
    global x, y, z
    x, y, z = symbols("x, y, z")

    while True:
        instruct_path = (
            os.path.dirname(os.path.abspath(__file__)) + "/texts/main_screen.txt"
        )
        main = open(instruct_path, mode="r")
        for line in main.readlines():
            line = line.rstrip()
            nprint(line)

        try:
            if date_option == "1":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d")
            elif date_option == "2":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")
        except:
            now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")

        nprint(f"\n(Time now is: {now})")
        nprint("(Current Screen: Main Screen)\n")
        cmd = input("Enter Command: ")

        if cmd == "1":
            derivacalc()

        elif cmd == "2":
            intecalc()

        elif cmd == "3":
            limcalc()

        elif cmd == "4":
            settings()

        elif cmd == "5":
            nprint("\nExiting Calc-ULTRA ... ... ...\n")
            break

        else:
            logging.warning(f'Invalid command:"{cmd}"\n')


"""
If you find this message, type 'hi' in the general discussions - sudoer-Huatao
"""


def derivacalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/derivacalc.txt"
    derivacalc = open(instruct_path, mode="r")
    for line in derivacalc.readlines():
        line = line.rstrip()
        nprint(line)

    while True:
        nprint("\n(Current Screen: DerivaCalc Main Screen)\n")
        cmd = input("Enter Command: ")

        if cmd == "1":
            nprint("\n(Current Screen: Derivative Screen)\n")
            global dfunction
            dfunction = input("Enter a function: ")
            order = input("Enter order of derivative calculation: ")
            derive(dfunction, order)

        elif cmd == "2":
            nprint("\n(Current Screen: Partial Derivative Screen)\n")
            function = input("Enter a function containing x and y or x and y and z: ")
            var = input("Enter variable to differentiate in respect to: ")
            if var != "x" and var != "y" and var != "z":
                logging.error("Variable to differentite in respect to is invalid.")
            else:
                order = input("Enter the order of partial derivative calculation: ")
                partial_derive(function, var, order)

        elif cmd == "3":
            nprint("\n(Current Screen: Implicit Derivative Screen)\n")
            circ = input(
                "Enter the left side of an equation containing x and y: (right side default as 0) "
            )
            order = input("Enter order of implicit derivative calculation: ")
            implicit_derive(circ, order)

        elif cmd == "4":
            nprint("\nExiting DerivaCalc ... ... ...")
            break

        else:
            logging.warning(f'Invalid command:"{cmd}"')


def intecalc():
    global dfunction

    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/intecalc.txt"
    intecalc = open(instruct_path, mode="r")
    for line in intecalc.readlines():
        line = line.rstrip()
        nprint(line)

    while True:
        nprint("(Current Screen: InteCalc Main Screen)\n")
        cmd = input("Enter Command: ")

        if cmd == "1":
            nprint("\n(Current Screen: Antiderivative Screen)\n")
            dfunction = input("Enter a function: ")
            antiderive(dfunction)

        elif cmd == "2":
            nprint("\n(Current Screen: Definite Integral Screen)\n")
            dfunction = input("Enter a function: ")
            lower_bound = input("\nEnter the lower bound: ")
            upper_bound = input("Enter the upper bound: ")
            print(definite_integrate(dfunction, lower_bound, upper_bound))

        elif cmd == "3":
            nprint("\n(Current Screen: Improper Integral Screen)\n")
            function = input("Enter a function: ")
            lower_bound = input("\nEnter the lower bound: ")
            upper_bound = input("Enter the upper bound: ")
            improper_integrate(function, lower_bound, upper_bound)

        elif cmd == "4":
            nprint("\nExiting InteCalc ... ... ...")
            break

        else:
            logging.warning(f'Invalid command: "{cmd}"')


def limcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/limcalc.txt"
    limcalc = open(instruct_path, mode="r")
    for line in limcalc.readlines():
        line = line.rstrip()
        nprint(line)

    while True:
        try:
            nprint("\n(Current Screen: LimCalc Main Screen)\n")
            cmd = input("Enter Command: ")

            if cmd == "1":
                nprint("\n(Current screen: Limit Screen)\n")
                expr = input("Enter an expression: ")
                value = input("Enter point of evaluation: ")
                normal_limit(expr, value)

            elif cmd == "2":
                nprint("\n(Current screen: One-sided Limit Screen)\n")
                expr = input("Enter an expression: ")
                value = input("Enter point of evaluation: ")
                direction = input("Enter direction of limit ('left' or 'right'): ")
                one_side_limit(expr, value, direction)

            elif cmd == "3":
                nprint("\nExiting LimCalc ... ... ...")
                break

            else:
                logging.warning(f'Invalid command: "{cmd}"')

        except:
            logging.error("UnknownError: An unknown error occured.")


def settings():
    settings_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/settings.txt"
    settings = open(settings_path, mode="r")

    for line in settings.readlines():
        line = line.rstrip()
        nprint(line)

    while True:
        nprint("\n(Current Screen: Settings Screen)\n")
        cmd = input("Enter Command: ")

        if cmd == "print":
            nprint("\n(Current Screen: Print Settings Screen)\n")

            global print_option
            print_option = input(
                'Set print mode: "p" (Sympy Pretty Print) or "n" (Normal Print): '
            )
            nprint(f'\nPrinting mode set to: "{print_option}"')

        elif cmd == "graph":
            nprint("\n(Current Screen: Graph Settings Screen)\n")

            global graph_option
            graph_option = input(
                'Set graph mode: "f" (Fixed graph view) or "a" (Adjusted graph view): '
            )
            nprint(f'\nGraph mode set to: "{graph_option}"')

        elif cmd == "date":
            nprint("\n(Current Screen: Date Settings Screen)\n")

            global date_option
            date_option = input(
                'Set date mode: "1" (YY/MM/DD) or "2" (YY/MM/DD/HH/MM/SS): '
            )
            nprint(f'\nDate mode set to: "{date_option}"')

        elif cmd == "exit":
            nprint("\nExiting settings ... ... ...")
            break

        else:
            logging.warning(f'Invalid command:"{cmd}"')


main()
