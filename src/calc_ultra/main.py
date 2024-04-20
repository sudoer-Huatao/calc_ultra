from sympy.core.numbers import pi, E, oo
from math import floor, ceil
import math as mt
from scipy.special import gamma, polygamma, erf
from rich import print
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    MofNCompleteColumn,
    TimeElapsedColumn,
    TextColumn,
)
from numpy import (
    linspace,
    exp,
    log,
    sqrt,
    abs,
    sin,
    cos,
    tan,
    arcsin,
    arccos,
    arctan,
    sinh,
    cosh,
    tanh,
    arcsinh,
    arccosh,
    arctanh,
    array,
    cross,
)
from numpy.linalg import norm, det
import numpy as np

# Sympy uses symbolic Pi and e, which cannot be graphed by matplotlib
# so later np.pi np.e are explicitly used instead.
from sympy import (
    diff,
    idiff,
    Integral,
    integrate,
    limit,
    solve,
    linsolve,
    pprint,
    simplify,
    symbols,
)
import matplotlib.pyplot as plt
import datetime, logging, random, os, time, warnings


# Disable auto Python warnings

warnings.filterwarnings("ignore")


def simp():
    print(
        '\n[bold bright_green](Current Screen: Simple Calculation Screen)[/bold bright_green]\n("q" to quit)\nEnter any expression to start:\n'
    )

    while True:
        expr = input()

        try:
            if expr != "q":
                # Easy to exit unlike VIM ;)
                result = eval(trig_rep(expr))
                print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                pprint(result)
                print()

            else:
                break

        except:
            print()
            logging.error(f'Could not parse: "{expr}"\n')


def derive(function: str, order: str):
    calc = replace_expr(function)

    if check_order(order) is True:
        global df
        df = diff(calc, x, order)

        print(
            f"\nDerivative of [bright_magenta]{trig_rep(function)}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        check_simp(df)

        df = trig_rep(str(df))

        print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
        show = input("(Exit the graph window when you are finished to continue) ")

        if show == "y":
            try:
                print()
                with Progress(
                    SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                    TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                    BarColumn(),
                    TimeElapsedColumn(),
                    transient=False,
                ) as progress:
                    task = progress.add_task("", total=100)

                    while not progress.finished:
                        progress.update(task, advance=2)
                        time.sleep(random.randint(2, 5) / 1000)

                x_array = linspace(-50, 50, 200000)

                def f(x):
                    return eval(trig_rep(calc))

                def dif(x):
                    return eval(df)

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

                print("\nExited graph.\n")

            except:
                plt.close()
                print("\n")
                logging.warning("Could not graph function.")
                print("\nExited graph.")


def partial_derive(function: str, var: str, order: str):
    calc = replace_expr(function)

    if check_order(order) is True:
        df = diff(calc, var, order)

        print(
            f"\nPartial derivative of [bright_magenta]{trig_rep(function)}[/bright_magenta] in respect to {var} of order {order} is:\n"
        )
        print_expr(df)

        check_simp(df)


def implicit_derive(circ: str, order: str):
    calc = replace_expr(circ)
    left = eval(calc[: calc.find("=")])
    right = eval(calc[calc.find("=") + 1 :])

    if check_order(order) is True:
        df = idiff(left - right, y, x, order)

        print(
            f"\nDerivative of [bright_magenta]{circ}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        if str(simplify(df, evaluate=False)) != str(df):
            print("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))


def antiderive(function: str):
    calc = replace_expr(function)
    F = Integral(calc, x).doit()

    if "Integral" in str(F):
        logging.warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nAntiderivative of [bright_magenta]{trig_rep(function)}[/bright_magenta] is:\n"
    )
    print_expr(F)

    check_simp(F)

    F = trig_rep(str(F))

    print("\n[bold]Don't forget to add a constant![/bold]\n")

    print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            print()
            with Progress(
                SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                BarColumn(),
                TimeElapsedColumn(),
                transient=False,
            ) as progress:
                task = progress.add_task("", total=100)

                while not progress.finished:
                    progress.update(task, advance=2)
                    time.sleep(random.randint(2, 5) / 1000)

            x_array = linspace(-100, 100, 200000)

            def f(x):
                return eval(trig_rep(calc))

            def af(x):
                return eval(F)

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

            print("\nExited graph.\n")

        except:
            plt.close()
            print("\n")
            logging.warning("Could not graph function.")
            print("\nExited graph.\n")


def def_int(function: str, low: str, up: str):
    x = symbols("x")
    calc = replace_expr(function)

    check_bound(low)

    clow = eval(replace_expr(low))

    check_bound(up)

    cup = eval(replace_expr(up))

    result = integrate(calc, (x, clow, cup))

    gup = eval(replace_bound(str(cup)))
    glow = eval(replace_bound(str(clow)))

    num_result = integrate(calc, (x, glow, gup)).evalf()
    # Composite functions usually do not have primitive antiderivatives
    # so calc-ultra is equipped with both symbolic and numerical answers.

    if (
        (str(result) == "nan")
        or ("I" in str(result))
        and ("Integral" not in str(result))
    ):
        logging.warning("Cannot compute integral because integral does not converge.")
        return ""

    if "Integral" not in str(result):
        print(
            f"\nCalculated integral of [bright_magenta]{trig_rep(function)}[/bright_magenta] from {low} to {up}. Final area is:\n"
        )
        print_expr(result)

    else:
        print("\nCannot express result symbolically.")

    nprint("\nNumeral evaluation/approximation:\n")
    print_expr(num_result)

    print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")

    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            print()
            with Progress(
                SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
                TextColumn("[bright_yellow]Loading graph...[/bright_yellow]"),
                BarColumn(),
                TimeElapsedColumn(),
                transient=False,
            ) as progress:
                task = progress.add_task("", total=100)

                while not progress.finished:
                    progress.update(task, advance=2)
                    time.sleep(random.randint(2, 5) / 1000)

            x_array = linspace(-100, 100, 200000)

            if "log" in function:
                x_array = linspace(
                    0.00000001,
                    100,
                    200000,
                )

            if "factorial" in function:
                x_array = linspace(
                    0,
                    100,
                    200000,
                )

            # This is just the gamma function shifted one unit left
            # Of course, the factorial is only defined for x >= 0

            def f(x):
                return eval(trig_rep(calc))

            title = "Shaded area beneath function"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_array, f(x_array), color="red", label="Function")

            plt.fill_between(
                x_array,
                f(x_array),
                where=[(x_array > clow) and (x_array < cup) for x_array in x_array],
                color="blue",
            )

            try:
                if graph_option == "f":
                    plt.axis([-7.5, 7.5, -7.5, 7.5])

                elif graph_option == "a":
                    # Adjusted graph view is sometimes better for
                    # large graphs with large bounds.
                    if (float(f(glow)) != 0) and (float(f(gup)) != 0):
                        plt.axis(
                            [
                                glow - 5,
                                gup + 5,
                                float(f(round(glow)))
                                - (float(f(round(glow))) + float(f(round(gup)))) / 2
                                - 1,
                                float(f(round(gup)))
                                + (float(f(round(glow))) + float(f(round(gup)))) / 2
                                + 1,
                            ]
                        )

                    elif (float(f(glow)) == 0) or (float(f(gup)) == 0):
                        plt.axis(
                            [glow - 5, gup + 5, -(gup - glow) / 2, (gup + glow) / 2]
                        )

            except:
                plt.axis([-7.5, 7.5, -7.5, 7.5])

            plt.legend(loc="lower left")
            plt.grid()
            plt.show()

            return "\nExited graph.\n"

        except:
            plt.close()
            print("\n")
            logging.warning("Could not graph function.")
            return "\nExited graph.\n"

    else:
        return "\n[bright_yellow]Exiting Definite Integral Screen ... ... ...[/bright_yellow]\n"


def improp_int(function: str, low: str, up: str):
    function = replace_expr(function)

    check_bound(low)

    clow = eval(replace_expr(low))

    check_bound(up)

    cup = eval(replace_expr(up))

    try:
        improper_area = Integral(
            function, (x, clow, cup)
        ).principal_value()  # Cauchy Principal Value

        if "Integral" not in str(improper_area):
            print(
                f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)

        else:
            print("Cannot compute improper integral.")

    except ValueError:
        improper_area = integrate(function, (x, clow, cup))

        if "Integral" not in str(improper_area):
            print(
                f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)

        else:
            print("Cannot compute improper integral.")

    print()


def double_int(function: str, out_low: str, out_up: str, in_low: str, in_up: str):
    function = replace_expr(function)

    check_bound(out_low)

    out_low = eval(replace_expr(out_low))

    check_bound(out_up)

    out_up = eval(replace_expr(out_up))

    in_low = eval(replace_expr(in_low))

    in_up = eval(replace_expr(in_up))

    out_up = eval(replace_bound(str(out_up)))
    out_low = eval(replace_bound(str(out_low)))
    in_up = eval(replace_bound(str(in_up)))
    in_low = eval(replace_bound(str(in_low)))

    result = integrate(function, (y, in_low, in_up), (x, out_low, out_up))

    if "Integral" in str(result):
        logging.warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nDouble integral of [bright_magenta]{function}[/bright_magenta] with inner bounds [bright_cyan]{in_low}[/bright_cyan] and [bright_cyan]{in_up}[/bright_cyan] and outer bounds {out_low} and {out_up} is:\n"
    )
    print_expr(result)
    print("")

    check_simp(result)


def lim(expr: str, value: str):
    expr = replace_expr(expr)

    check_bound(value)

    value = float(eval(replace_bound(value)))

    l = limit(expr, x, value)

    if "Limit" in str(l):
        logging.warning("Cannot compute limit.")
        return ""

    if limit(expr, x, value, "+") != limit(expr, x, value, "-"):
        print(
            "\nThe limit does not exist (i.e., the limit approaching from the right does not equal the limit approaching from the left)."
        )
        print(
            "Use the one-side limit calculation of LimCalc to calculate the limit at one side.\n"
        )

    else:
        print(
            f"\nLimit of [bright_magenta]{expr}[/bright_magenta] as x approaches {value} is:\n"
        )
        print_expr(l)
        check_simp(l)


def side_lim(expr: str, value: str, direction: str):
    expr = replace_expr(expr)

    check_bound(value)

    value = float(eval(replace_bound(value)))

    if direction == "left" or direction == "Left":
        direction_sign = "-"

    elif direction == "right" or direction == "Right":
        direction_sign = "+"

    else:
        print()
        logging.error("\nTypeError: Direction is neither right or left.")
        return ""

    l = limit(expr, x, value, dir=direction_sign)

    if "Limit" in str(l):
        logging.warning("\nCannot compute limit.")
        return ""

    print(
        f"\nLimit of [bright_magenta]{expr}[/bright_magenta] as x approaches {value} from the [bright_cyan]{direction}[/bright_cyan] is:\n"
    )
    print_expr(l)


def eq_solve(mode: int):
    eq_list = []

    if mode == 1:
        eq1 = input("\nEnter equation: ")
        left = eval(eq1[: eq1.find("=")])
        right = eval(eq1[eq1.find("=") + 1 :])
        eq_set = solve(left - right)
        if len(eq_set) == 0:
            print()
            logging.error("UnknownError: Cannot solve equation")
        else:
            print("\nx:\n")
            for i in range(0, len(eq_set)):
                pprint(eq_set[i])
                print()

    elif mode == 2:
        eq1 = input("Enter first equation: ")
        left1 = eval(eq1[: eq1.find("=")])
        right1 = eval(eq1[eq1.find("=") + 1 :])

        eq2 = input("Enter second equation: ")
        left2 = eval(eq2[: eq2.find("=")])
        right2 = eval(eq2[eq2.find("=") + 1 :])

        eq_set = str(linsolve((left1 - right1, left2 - right2), (x, y), (-1, 1)))
        eqs = eq_set.strip("{").strip("}").strip("(").strip(")")
        for value in eqs:
            eq_list.append(value)

        result = "".join(eq_list)
        print("\nx, y:\n")
        pprint(result)

    elif mode == 3:
        eq1 = input("Enter equation 1: ")
        left1 = eval(eq1[: eq1.find("=")])
        right1 = eval(eq1[eq1.find("=") + 1 :])

        eq2 = input("Enter equation 2: ")
        left2 = eval(eq2[: eq2.find("=")])
        right2 = eval(eq2[eq2.find("=") + 1 :])

        eq3 = input("Enter equation 3: ")
        left3 = eval(eq3[: eq3.find("=")])
        right3 = eval(eq3[eq3.find("=") + 1 :])

        eq_set = str(
            linsolve(
                (left1 - right1, left2 - right2, left3 - right3), (x, y, z), (-1, 1)
            )
        )
        eqs = eq_set.strip("{").strip("}").strip("(").strip(")")
        for value in eqs:
            eq_list.append(value)

        result = "".join(eq_list)
        print("\nx, y, z:\n")
        pprint(result)


def check_simp(expr) -> bool:
    if str(simplify(expr, evaluate=False)) != str(expr):
        print("\nSimplify/rewrite:\n")
        print_expr(simplify(expr, evaluate=False))

    else:
        return False


def check_order(order: str) -> bool:
    if ("." in order) or (order.isnumeric() == False):
        print()
        logging.error(
            "OrderError: Order of derivative calculation is not a valid number."
        )
        return False

    else:
        return True


def check_bound(bound: str):
    if (
        (bound.isnumeric() is False)
        and ("pi" not in bound)
        and ("e" not in bound)
        and ("E" not in bound)
        and ("-" not in bound)
        and ("." not in bound)
        and ("sqrt" not in bound)
        and ("oo" not in bound)
        and ("/" not in bound)
    ):
        print()
        logging.error("TypeError: Integration bound is not a number.")


def replace_expr(expr: str) -> str:
    expr = expr.strip(" ")

    if "^" in expr:
        expr = expr.replace("^", "**")

    if "e" in expr:
        expr = expr.replace("e", "E")

    if "E**" in expr:
        expr = expr.replace("E**", "exp")

    if "ln" in expr:
        expr = expr.replace("ln", "log")

    if "arc" in expr:
        expr = expr.replace("arc", "a")

    if "abs" in expr:
        expr = expr.replace("abs", "Abs")

    return expr


def replace_bound(bound: str) -> str:
    if "pi" in bound:
        bound = bound.replace("pi", str(np.pi))

    if "E" in bound:
        bound = bound.replace("E", str(np.e))

    return bound


def factorial(x):
    return gamma(x + 1)


def trig_rep(function: str) -> str:
    # Sympy and Numpy trig functions are vastly different
    # and uses different prefixes. Thus this replacement
    # algorithm is needed.
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

    if "csc" in function:
        function = function.replace("csc", "1/sin")

    # Unfortunately, Numpy does not have an implemented csc,
    # sec, cot, csch, etc. so we have to implement our own.

    if "sec" in function:
        function = function.replace("sec", "1/cos")

    if "cot" in function:
        function = function.replace("cot", "1/tan")

    if "csch" in function:
        function = function.replace("csch", "1/sinh")

    if "sech" in function:
        function = function.replace("sech", "1/cosh")

    if "coth" in function:
        function = function.replace("coth", "1/tanh")

    if "Abs" in function:
        function = function.replace("Abs", "fabs")

    if "log" in function and ",x)" in function:
        function = function.replace("log", "mt.log")

    if "x" not in function:
        function = "0 * x + " + function

    return function


def print_expr(text: str):
    printing_methods = {"p": lambda t: pprint(text), "n": lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods["p"](text)


def nprint(text: str):
    print(text)
    time.sleep(0.04)


def main():
    with Progress(
        SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
        TextColumn("[green]Handling imports[/green]..."),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        transient=False,
    ) as progress:
        task = progress.add_task("", total=50)

        while not progress.finished:
            progress.update(task, advance=1)
            time.sleep(random.randint(2, 5) / 100)

    global x, y, z
    x, y, z = symbols("x, y, z")

    while True:
        instruct_path = (
            os.path.dirname(os.path.abspath(__file__))
            + "/texts/main_screen.txt"
            # TODO: make the PATH compatible with Windows
        )
        main = open(instruct_path, mode="r")
        for line in main.readlines():
            line = line.rstrip()
            if ("---" in line) or ("|" in line):
                nprint("[gold1]" + line + " [/gold1]")
            else:
                nprint(line)

        try:
            if date_option == "1":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d")
            elif date_option == "2":
                now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")
        except:
            now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")

        print(f"\n(Time now is: {now})\n")
        print("[bold bright_green](Current Screen: Main Screen)[/bold bright_green]\n")
        print("[bright_magenta]Enter Command: [/bright_magenta]", end="")
        cmd = input()

        if cmd == "1":
            simp()

        elif cmd == "2":
            derivacalc()

        elif cmd == "3":
            intecalc()

        elif cmd == "4":
            limcalc()

        elif cmd == "5":
            algcalc()

        elif cmd == "6":
            settings()

        elif cmd == "7":
            nprint("\n[bright_yellow]Exiting Calc-ULTRA ... ... ...[/bright_yellow]\n")
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
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: DerivaCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_magenta]Enter Command: [/bright_magenta]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current Screen: Derivative Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                order = input("Enter order of derivative calculation: ")
                derive(function, order)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Partial Derivative Screen)[/bold bright_green]\n"
                )
                function = input(
                    "Enter a function containing x and y or x and y and z: "
                )
                var = input("Enter variable to differentiate in respect to: ")
                if var != "x" and var != "y" and var != "z":
                    print()
                    logging.error("Variable to differentite in respect to is invalid.")
                else:
                    order = input("Enter the order of partial derivative calculation: ")
                    partial_derive(function, var, order)

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current Screen: Implicit Derivative Screen)[/bold bright_green]\n"
                )
                circ = input("Enter an equation containing x and y:")
                order = input("Enter order of implicit derivative calculation: ")
                implicit_derive(circ, order)

            elif cmd == "4":
                print("\n[bright_yellow]Exiting DerivaCalc ... ... ...[/bright_yellow]")
                break

            else:
                logging.warning(f'Invalid command:"{cmd}"')
        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def intecalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/intecalc.txt"
    intecalc = open(instruct_path, mode="r")
    for line in intecalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "[bold bright_green](Current Screen: InteCalc Main Screen)[/bold bright_green]\n"
            )
            print("[bright_magenta]Enter Command: [/bright_magenta]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current Screen: Antiderivative Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                antiderive(function)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Definite Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                lower_bound = input("\nEnter the lower bound: ")
                upper_bound = input("Enter the upper bound: ")
                print(def_int(function, lower_bound, upper_bound))

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current Screen: Improper Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                lower_bound = input("\nEnter the lower bound: ")
                upper_bound = input("Enter the upper bound: ")
                improp_int(function, lower_bound, upper_bound)

            elif cmd == "4":
                nprint(
                    "\n[bold bright_green](Current Screen: Double Integral Screen)[/bold bright_green]\n"
                )
                function = input("Enter a function: ")
                outer_low = input("\nEnter the lower outer bound: ")
                outer_up = input("Enter the upper outer bound: ")
                inner_low = input("\nEnter the lower inner bound: ")
                inner_up = input("Enter the upper inner bound: ")
                double_int(function, outer_low, outer_up, inner_low, inner_up)

            elif cmd == "5":
                print("\n[bright_yellow]Exiting InteCalc ... ... ...[/bright_yellow]")
                break

            else:
                logging.warning(f'Invalid command: "{cmd}"')
        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def limcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/limcalc.txt"
    limcalc = open(instruct_path, mode="r")
    for line in limcalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: LimCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current screen: Limit Screen)[/bold bright_green]\n"
                )
                expr = input("Enter an expression: ")
                value = input("Enter point of evaluation: ")
                lim(expr, value)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current screen: One-sided Limit Screen)[/bold bright_green]\n"
                )
                expr = input("Enter an expression: ")
                value = input("Enter point of evaluation: ")
                direction = input("Enter direction of limit ('left' or 'right'): ")
                side_lim(expr, value, direction)

            elif cmd == "3":
                print("\n[bright_yellow]Exiting LimCalc ... ... ...[/bright_yellow]")
                break

            else:
                logging.warning(f'Invalid command: "{cmd}"')

        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def algcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + "/texts/algcalc.txt"
    algcalc = open(instruct_path, mode="r")
    for line in algcalc.readlines():
        line = line.rstrip()
        if ("---" in line) or ("|" in line):
            nprint("[gold1]" + line + " [/gold1]")
        else:
            nprint(line)
    print()

    while True:
        try:
            nprint(
                "\n[bold bright_green](Current Screen: AlgCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                nprint(
                    "\n[bold bright_green](Current screen: Equation Solver Screen)[/bold bright_green]\n"
                )
                print(
                    "Enter mode: 1 for one set equation, 2 for two, and 3 for three: ",
                    end="",
                )
                mode = input()
                eq_solve(mode)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current screen: Vector Calculation Screen)[/bold bright_green]\n"
                )

                print("Write a vector like [1, 2, 3], then perform operations!")
                print(
                    "\n- Use [bright_magenta]@[/bright_magenta] to calculate the dot product"
                )
                print("- cross(smth, smth) to calculate the cross product")
                print(
                    "- A pair of [bright_magenta]< >[/bright_magenta]s encasing a vector to calculate it's norm!"
                )
                print('Enter any expression to start ("q" to quit):\n')

                while True:
                    expr = input()
                    try:
                        if expr != "q":
                            expr = (
                                (
                                    (
                                        (expr.replace("[", "array([")).replace(
                                            "]", "])"
                                        )
                                    ).replace("<", "norm(")
                                ).replace(">", ")")
                            ).strip(" ")
                            result = eval(trig_rep(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        logging.error(f'Could not parse: "{expr}"\n')

                print()

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current screen: Matrix Calculation Screen)[/bold bright_green]\n"
                )

                print("Write a matrix like [[1, 2], [3, 4]], then perform operations!")
                print(
                    "- Use [bright_magenta]@[/bright_magenta] to calculate the dot product"
                )
                print(
                    "- A pair of [bright_magenta]| |[/bright_magenta] to calculate the determinant."
                )
                print('Enter any expression to start ("q" to quit):\n')

                while True:
                    expr = input()
                    try:
                        if expr != "q":
                            expr = (
                                (
                                    (
                                        (expr.replace("[[", "array([[")).replace(
                                            "]]", "]])"
                                        )
                                    ).replace("|a", "det(a")
                                ).replace(")|", "))")
                            ).strip(" ")
                            result = eval(trig_rep(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]\n")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        print("\n")
                        logging.error(f'Could not parse: "{expr}"\n')

            elif cmd == "4":
                print("\n[bright_yellow]Exiting AlgCalc ... ... ...[/bright_yellow]")
                break

            else:
                logging.warning(f'Invalid command: "{cmd}"')

        except:
            print("\n")
            logging.error("UnknownError: An unknown error occured.")
            nprint("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def settings():
    while True:
        settings_path = (
            os.path.dirname(os.path.abspath(__file__)) + "/texts/settings.txt"
        )
        settings = open(settings_path, mode="r")

        for line in settings.readlines():
            line = line.rstrip()
            if ("---" in line) or ("|" in line):
                nprint("[gold1]" + line + " [/gold1]")
            else:
                nprint(line)
        print("\n[bold green](Current Screen: Settings Screen)[/bold green]\n")
        print("[bright_magenta]Enter Command: [/bright_magenta]", end="")
        cmd = input()

        if cmd == "print":
            print(
                "\n[bold bright_green](Current Screen: Print Settings Screen)[/bold bright_green]\n"
            )

            global print_option
            print_option = input(
                'Set print mode: "p" (Sympy Pretty Print) or "n" (Normal Print): '
            )
            print(f'\nPrinting mode set to: "{print_option}"')

        elif cmd == "graph":
            print(
                "\n[bold bright_green](Current Screen: Graph Settings Screen)[/bold bright_green]\n"
            )

            global graph_option
            graph_option = input(
                'Set graph mode: "f" (Fixed graph view) or "a" (Adjusted graph view): '
            )
            print(f'\nGraph mode set to: "{graph_option}"')

        elif cmd == "date":
            print(
                "\n[bold bright_green](Current Screen: Date Settings Screen)[/bold bright_green]\n"
            )

            global date_option
            date_option = input(
                'Set date mode: "1" (YY/MM/DD) or "2" (YY/MM/DD/HH/MM/SS): '
            )
            print(f'\nDate mode set to: "{date_option}"')

        elif cmd == "style":
            print(
                "\n[bold bright_green](Current Screen: Graph Style Settings Screen)[/bold bright_green]\n"
            )

            print("You currently have these available styles:\n")

            style_list = plt.style.available
            styles = ", ".join(style_list)
            nprint("[bright_yellow]" + styles + " [/bright_yellow]")

            style = input('\nChoose a style to apply (type "default" to reset): ')

            plt.style.use(style)

            print(f'Graph style set to "{style}".')

        elif cmd == "exit":
            print("\n[bright_yellow]Exiting settings ... ... ...[/bright_yellow]")
            break

        else:
            logging.warning(f'Invalid command:"{cmd}"')


main()
