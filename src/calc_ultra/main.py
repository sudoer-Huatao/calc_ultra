from sympy.core.numbers import pi, E, oo
from math import floor, ceil
from scipy.special import gamma, polygamma, erf
from rich import print
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
)
import numpy as np

# Sympy uses symbolic Pi and e, which cannot be graphed by matplotlib
# so later np.pi np.e are explicitly used instead.
from sympy import (
    diff,
    idiff,
    Integral,
    integrate,
    limit,
    pprint,
    simplify,
    symbols,
)
import matplotlib.pyplot as plt
import datetime, logging, os, time, warnings


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
                print("\nResult:")
                if "pi" in expr:
                    expr = expr.replace("pi", str(np.pi))

                result = eval(replace_bound(expr))
                pprint(result)
                print()

            else:
                break

        except:
            logging.error(f'Could not parse: "{expr}"\n')


def derive(function, order):
    function = replace_expr(function)

    if check_order(order) is True:
        global df
        df = diff(function, x, order)

        nprint(
            f"\nDerivative of [bright_magenta]{function}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        if check_simp(df) is True:
            nprint("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))

        df = trig_rep(str(df))

        print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
        show = input("(Exit the graph window when you are finished to continue) ")

        if show == "y":
            try:
                print(
                    "\n[bright_yellow]Loading graph. Might take some time on first startup ...[/bright_yellow]"
                )
                x_array = linspace(-50, 50, 200000)

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
                logging.warning("Could not graph function.")
                print("\nExited graph.")


def partial_derive(function, var, order):
    function = replace_expr(function)

    if check_order(order) is True:
        df = diff(function, var, order)

        print(
            f"\nPartial derivative of [bright_magenta]{function}[/bright_magenta] in respect to {var} of order {order} is:\n"
        )
        print_expr(df)

        if check_simp(df) is True:
            print("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))


def implicit_derive(circ, order):
    circ = replace_expr(circ)

    if check_order(order) is True:
        df = idiff(eval(circ), y, x, int(order))

        print(
            f"\nDerivative of [bright_magenta]{circ}[/bright_magenta] with order {order} is:\n"
        )
        print_expr(df)

        if str(simplify(df, evaluate=False)) != str(df):
            print("\nSimplify/rewrite:\n")
            print_expr(simplify(df, evaluate=False))


def antiderive(function):
    function = replace_expr(function)

    global F
    F = Integral(function, x).doit()

    if "Integral" in str(F):
        logging.warning("Cannot compute integral.\n")
        return ""

    print(f"\nAntiderivative of [bright_magenta]{function}[/bright_magenta] is:\n")
    print_expr(F)

    if check_simp(F) is True:
        print("\nSimplify/rewrite:\n")
        print_expr(simplify(F, evaluate=False))

    F = trig_rep(str(F))

    print("\n[bold]Don't forget to add a constant![/bold]\n")

    print("\n[bright_yellow]Show graph of area? (y/n)[/bright_yellow]")
    show = input("(Exit the graph window when you are finished to continue) ")

    if show == "y":
        try:
            print(
                "\n[bright_yellow]Loading graph. Might take some time on first startup ...[/bright_yellow]"
            )

            x_array = linspace(-100, 100, 200000)

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
            logging.warning("Could not graph function.")
            print("\nExited graph.\n")


def def_int(function, low, up):
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
            f"\nCalculated integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
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
            print(
                "\n[bright_yellow]Loading graph. Might take some time on first startup ...[/bright_yellow]"
            )

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
                    # Adjusted graph view is sometimes better for
                    # large graphs with large bounds.
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

            return "\nExited graph.\n"

        except:
            plt.close()
            logging.warning("Could not graph function.")
            return "\nExited graph.\n"

    else:
        return "\n[bright_yellow]Exiting Definite Integral Screen ... ... ...[/bright_yellow]\n"


def improp_int(function, low, up):
    function = replace_expr(function)

    if check_bound(low) is False:
        return ""

    low = eval(low)

    if check_bound(up) is False:
        return ""

    up = eval(up)

    try:
        improper_area = Integral(
            function, (x, low, up)
        ).principal_value()  # Cauchy Principal Value

        if "Integral" not in str(improper_area):
            print(
                f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
            )
            print_expr(improper_area)

        else:
            print("Cannot compute improper integral.")

        print()

    except ValueError:
        logging.warning("ValueError: Singularity while computing improper integral.\n")


def double_int(function, out_low, out_up, in_low, in_up):
    function = replace_expr(function)

    if check_bound(out_low) is False:
        return ""

    out_low = eval(out_low)

    if check_bound(out_up) is False:
        return ""

    out_up = eval(out_up)

    in_low = eval(in_low)

    in_up = eval(in_up)

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

    if check_simp(result) is True:
        print("\nSimplify/rewrite:\n")
        print_expr(simplify(result, evaluate=False))


def lim(expr, value):
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
        if check_simp(l) is True:
            print("\nSimplify/rewrite:\n")
            print_expr(simplify(l, evaluate=False))


def side_lim(expr, value, direction):
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

    print(
        f"\nLimit of [bright_magenta]{expr}[/bright_magenta] as x approaches {value} from the [bright_cyan]{direction}[/bright_cyan] is:\n"
    )
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
        logging.error("TypeError: Integration bound is not a number.")
        return False

    else:
        return True


def replace_expr(expr):
    if "^" in expr:
        expr = expr.replace("^", "**")

    if "arc" in expr:
        expr = expr.replace("arc", "a")
        print(
            f'[bold bright_red]Showing results for "{expr}". (Did you make a mistake?)[/bold bright_red]'
        )

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
    if "x" not in F:
        final = eval("0 * x + " + F)

    else:
        final = eval(F)

    return final


def factorial(x):
    return gamma(x + 1)


def trig_rep(function):
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

    if "pi" in function:
        function = function.replace("pi", str(np.pi))

    if "E" in function:
        function = function.replace("E", str(np.e))

    return function


def print_expr(text):
    printing_methods = {"p": lambda t: pprint(text), "n": lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods["p"](text)


def err():
    nprint(
        "\n[bold bright_red]Check if your input is valid. You might have made the following mistakes:[/bold bright_red]"
    )
    nprint(
        '[bold bright_red] - Inverse trigs: "arc" instead of "a". "asin(x)" is correct, "arcsin(x)" is not.[/bold bright_red]'
    )
    nprint('[bold bright_red] - Base of log: "e" instead of "E".[/bold bright_red]\n')


def nprint(text):
    print(text)
    time.sleep(0.04)


def main():
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
            derivacalc()

        elif cmd == "2":
            intecalc()

        elif cmd == "3":
            limcalc()

        elif cmd == "4":
            simp()

        elif cmd == "5":
            settings()

        elif cmd == "6":
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
                global dfunction
                dfunction = input("Enter a function: ")
                order = input("Enter order of derivative calculation: ")
                derive(dfunction, order)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Partial Derivative Screen)[/bold bright_green]\n"
                )
                function = input(
                    "Enter a function containing x and y or x and y and z: "
                )
                var = input("Enter variable to differentiate in respect to: ")
                if var != "x" and var != "y" and var != "z":
                    logging.error(
                        "[bold bright_red]Variable to differentite in respect to is invalid.[/bold bright_red]"
                    )
                else:
                    order = input("Enter the order of partial derivative calculation: ")
                    partial_derive(function, var, order)

            elif cmd == "3":
                nprint(
                    "\n[bold bright_green](Current Screen: Implicit Derivative Screen)[/bold bright_green]\n"
                )
                circ = input(
                    "Enter the left side of an equation containing x and y: (right side default as 0) "
                )
                order = input("Enter order of implicit derivative calculation: ")
                implicit_derive(circ, order)

            elif cmd == "4":
                print("\n[bright_yellow]Exiting DerivaCalc ... ... ...[/bright_yellow]")
                break

            else:
                logging.warning(f'Invalid command:"{cmd}"')
        except:
            logging.error("UnknownError: An unknown error occured.")
            err()


def intecalc():
    global dfunction

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
                dfunction = input("Enter a function: ")
                antiderive(dfunction)

            elif cmd == "2":
                nprint(
                    "\n[bold bright_green](Current Screen: Definite Integral Screen)[/bold bright_green]\n"
                )
                dfunction = input("Enter a function: ")
                lower_bound = input("\nEnter the lower bound: ")
                upper_bound = input("Enter the upper bound: ")
                print(def_int(dfunction, lower_bound, upper_bound))

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
            logging.error(
                "UnknownError: An unknown error occured. Try the following things:"
            )
            err()


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
            logging.error(
                "UnknownError: An unknown error occured. Try the following things:"
            )
            err()


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
