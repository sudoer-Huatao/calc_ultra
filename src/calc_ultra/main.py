###########
# Imports #
###########


from rich import print

print(
    "\n[bright_yellow]Initializing Calc-Ultra. Might take some time on first startup ...[/bright_yellow]"
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
    gamma,
    lowergamma,
    uppergamma,
    zeta,
)
from sympy.core.numbers import pi, E, I, oo
import scipy.special as ss
from scipy.special import polygamma, gammainc, gammaincc, erf, erfc
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    MofNCompleteColumn,
    TimeElapsedColumn,
    TextColumn,
)
import matplotlib.pyplot as plt
from logging import warning, error
from math import floor, ceil
from os.path import dirname, abspath, exists
from prompt_toolkit import prompt
from prompt_toolkit.key_binding import KeyBindings
from random import randint
from readline import parse_and_bind
from time import sleep
from warnings import filterwarnings
import datetime


###########
# Presets #
###########

filterwarnings("ignore")

parse_and_bind("tab: complete")

x, y, z, t, m = symbols("x, y, z, t, m")

history = [""]  # Stores calculation history

key_bindings = KeyBindings()

idx = 1

save_path = dirname(abspath(__file__)) + "/fig.pdf"

expr_replace = [
    ("E", "E1"),
    ("e", "E1"),
    # The extra 1 at the end is to recognize whether the "e"
    # is part of an expression or a constant on it's own.
    ("i", "I"),
    ("^", "**"),
    ("E1**", "exp"),
    ("E1xp", "exp"),
    ("E1r", "er"),
    ("sIn", "sin"),
    ("pI", "pi"),
    ("Ia", "ia"),
    ("ln", "log"),
    ("arc", "a"),
    ("pi", "3.141592653589793"),
    ("E1", "2.718281828459045"),
]

graph_replace = [
    ("asin", "arcsin"),
    ("acos", "arccos"),
    ("atan", "arctan"),
    ("csc", "1/sin"),
    ("sec", "1/cos"),
    ("cot", "1/tan"),
    ("gamma", "ss.gamma"),
    ("polyss.gamma", "polygamma"),
    ("lowergamma", "gammainc"),
    ("uppergamma", "gammaincc"),
]


# Custom up/down key bindings


@key_bindings.add("up")
def _(event):
    global idx

    try:
        if idx != -1:
            idx -= 1
            buffer = event.app.current_buffer
            buffer.text = history[idx]
            buffer.cursor_position = len(buffer.text)
            idx = idx % len(history)
        else:
            idx = 0
            buffer = event.app.current_buffer
            buffer.text = history[idx]
            buffer.cursor_position = len(buffer.text)

    except:
        pass


@key_bindings.add("down")
def _(event):
    global idx

    try:
        if idx != len(history):
            idx += 1
            idx = idx % len(history)
            buffer = event.app.current_buffer
            buffer.text = history[idx]
            buffer.cursor_position = len(buffer.text)

        else:
            idx = len(history) - 1
            buffer = event.app.current_buffer
            buffer.text = history[idx]
            buffer.cursor_position = len(buffer.text)

    except:
        pass


# Control + C to clear history


@key_bindings.add("c-c")
def _(event):
    global history
    history = [""]


#####################
# Calculation funcs #
#####################


def simp():
    print(
        '\n[bold bright_green](Current Screen: Simple Calculation Screen)[/bold bright_green]\n\n("q" to quit)\n\nEnter any expression:\n'
    )

    while True:
        expr = prompt(key_bindings=key_bindings)

        try:
            if expr != "q":
                # Easy to exit unlike VIM ;)
                expr = expr.replace(" ", "")

                if "?=" in expr:
                    left = replace_expr(expr[: expr.find("?=")])
                    right = replace_expr(expr[expr.find("?=") + 2 :])

                    if "x" in left or "x" in right:
                        equals = "True."

                        left = str(left)
                        right = str(right)

                        for i in range(-100, 100):
                            left_result = simplify(
                                left.replace("x", str(i)).replace("e" + str(i), "ex")
                            ).evalf()
                            right_result = simplify(right.replace("x", str(i))).evalf()

                            if (
                                "I" not in str(left_result)
                                and left_result < 0.000001
                                and left_result > -0.000001
                            ):
                                left_result = 0

                            if (
                                "I" not in str(right_result)
                                and right_result < 0.000001
                                and right_result > -0.000001
                            ):
                                right_result = 0

                            if str(left_result).rstrip("*I") != str(
                                right_result
                            ).rstrip("*I"):
                                equals = "False."
                                print(left_result, right_result)
                                break

                        print(f"\n{equals}\n")

                    else:
                        left = replace_calc(left)
                        right = replace_calc(right)

                        if str(left).endswith("0") and left != 0:
                            left = str(left).rstrip("0")

                        if str(right).endswith("0") and right != 0:
                            right = str(right).rstrip("0")

                        if left == right:
                            print("\nTrue.\n")

                        else:
                            print(
                                "\nFalse.\n\nLeft side: "
                                + left
                                + ". Right side: "
                                + right
                                + ".\n"
                            )

                else:
                    result = replace_expr(expr)

                    result = replace_calc(result)

                    if result == None:
                        error("Could not complete calculation.")

                    else:
                        if str(result).endswith("0") and result != 0:
                            result = eval(str(result).rstrip("0"))

                        if history[len(history) - 1] != str(result) and result != None:
                            history.append(str(result))

                        print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                        pprint(result)
                        print()

            else:
                break

        except:
            print()
            error(f'Could not parse: "{expr}".\n')


def derive(function: str, order: str):
    r"""Calculates the derivative of a function.

    Explanation
    ===========

    Takes a function and an order of differentiation and
    returns the derivative of the function with the given order.
    Graph option available.
    """

    calc = replace_expr(function)

    check_order(order)

    global df
    df = diff(calc, x, order)

    print(
        f"\nDerivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] with order {order} is:\n"
    )
    print_expr(df)

    df = replace_expr(str(df))
    check_simp(df)

    if history[len(history) - 1] != df:
        history.append(df)

    df = replace_graph(df)

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
                    sleep(randint(2, 5) / 1000)

            x_arr = linspace(-100, 100, 200000)

            def f(x: array):
                return eval(replace_graph(calc))

            def dif(x: array):
                return eval(df)

            title = "Function (red) and derivative (blue)"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_arr, f(x_arr), color="red", label="Function")
            plt.plot(x_arr, dif(x_arr), color="blue", label="Derivative")
            plt.axis([-7.5, 7.5, -7.5, 7.5])
            plt.legend(loc="lower left")
            plt.grid()

            download = input("Download graph? (y/n) ")

            if download == "y":
                plt.savefig(save_path)
                print("\nDownloaded graph.")

            plt.show()
            print("\nExited graph.")

        except:
            plt.close()
            warning("Could not graph function.")
            print("\nExited graph.")


def partial_derive(function: str, var: str, order: str):
    r"""Calculates the partial derivative of a function.

    Explanation
    ===========

    Takes a function, a variable, and an order of differentiation and
    returns the partial derivative of the function in respect
    to the variable with the order.
    """

    calc = replace_expr(function)

    check_order(order)

    df = diff(calc, var, order)

    print(
        f"\nPartial derivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] in respect to {var} of order {order} is:\n"
    )
    print_expr(df)

    df = replace_expr(str(df))
    check_simp(df)

    if history[len(history) - 1] != df:
        history.append(df)


def implicit_derive(circ: str, order: str):
    r"""Calculates the implicit derivative of an equation.

    Explanation
    ===========

    Takes an equation and an order of differentiation and
    returns the implicit derivative of an equation with the given order.
    """

    calc = replace_expr(circ)
    left = eval(calc[: calc.find("=")])
    right = eval(calc[calc.find("=") + 1 :])

    check_order(order)

    df = idiff(left - right, y, x, order)

    print(
        f"\nDerivative of [bright_magenta]{circ}[/bright_magenta] with order {order} is:\n"
    )
    print_expr(df)

    df = replace_expr(str(df))
    check_simp(df)

    if history[len(history) - 1] != df:
        history.append(df)


def antiderive(function: str):
    r"""Calculates the antiderivative of a function.

    Explanation
    ===========

    Takes a function and returns the antiderivative of the function.
    Graph option available.
    """

    calc = replace_expr(function)
    F = Integral(calc, x).doit()

    if "Integral" in str(F):
        warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nAntiderivative of [bright_magenta]{replace_graph(function)}[/bright_magenta] is:\n"
    )
    print_expr(F)

    F = replace_expr(str(F))
    check_simp(F)

    if history[len(history) - 1] != F:
        history.append(F)

    F = replace_graph(F)

    print("\n[bold]Don't forget to add a constant![/bold]")

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
                    sleep(randint(2, 5) / 1000)

            x_arr = linspace(-100, 100, 200000)

            def f(x: array):
                return eval(replace_graph(calc))

            def af(x: array):
                return eval(F)

            title = "Function (red) and antiderivative (blue, C = 0)"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_arr, f(x_arr), color="red", label="Function")
            plt.plot(x_arr, af(x_arr), color="blue", label="Antiderivative")
            plt.axis([-7.5, 7.5, -7.5, 7.5])
            plt.legend(loc="lower left")
            plt.grid()
            download = input("Download graph? (y/n) ")

            if download == "y":
                plt.savefig(save_path)
                print("\nDownloaded graph.")

            plt.show()

            print("\nExited graph.")

        except:
            plt.close()
            warning("Could not graph function.")
            print("\nExited graph.")


def def_int(function: str, low: str, up: str):
    r"""Calculates the definite integral of a function.

    Explanation
    ===========

    Takes a function and a lower and upper bound and
    returns the definite integral from the lower bound to
    the upper bound of the function.
    Graph option available.
    """

    calc = replace_expr(function)

    check_num(low)

    calc_low = eval(replace_expr(low))

    check_num(up)

    calc_up = eval(replace_expr(up))

    result = integrate(calc, (x, calc_low, calc_up))

    graph_up = eval(replace_expr(str(calc_up)))
    graph_low = eval(replace_expr(str(calc_low)))

    num_result = integrate(calc, (x, graph_low, graph_up)).evalf()
    # Composite functions usually do not have primitive antiderivatives
    # so calc-ultra is equipped with both symbolic and numerical answers.

    if (
        (str(result) == "nan")
        or ("i" in str(result))
        and ("Integral" not in str(result))
    ):
        warning("Cannot compute integral because integral does not converge.")
        return ""

    if "Integral" not in str(result):
        print(
            f"\nCalculated integral of [bright_magenta]{replace_graph(function)}[/bright_magenta] from {low} to {up}. Final area is:\n"
        )
        print_expr(result)
        result = replace_expr(str(result))

        if history[len(history) - 1] != result:
            history.append(result)

    else:
        print("\nCannot express result symbolically.")

    print("\nNumeral evaluation/approximation:\n")
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
                    sleep(randint(2, 5) / 1000)

            x_arr = linspace(-100, 100, 200000)

            if "log" in function:
                x_arr = linspace(
                    0.0001,
                    100,
                    200000,
                )

            # The factorial is just the gamma function shifted one unit left
            # Of course, the factorial is only defined for x >= 0

            if "factorial" in function:
                x_arr = linspace(
                    0,
                    100,
                    200000,
                )

            def f(x: array):
                return eval(replace_graph(calc))

            title = "Shaded area beneath function"
            plt.title(title)
            plt.xlabel("x", weight="bold")
            plt.ylabel("y", rotation=0, weight="bold")
            plt.plot(x_arr, f(x_arr), color="red", label="Function")

            plt.fill_between(
                x_arr,
                f(x_arr),
                where=[(x_arr > calc_low) and (x_arr < calc_up) for x_arr in x_arr],
                color="blue",
            )

            try:
                if graph_option == "f":
                    plt.axis([-7.5, 7.5, -7.5, 7.5])

                elif graph_option == "a":
                    # Adjusted graph view is sometimes better for
                    # large graphs with large bounds.
                    if (float(f(graph_low)) != 0) and (float(f(graph_up)) != 0):
                        plt.axis(
                            [
                                graph_low - 5,
                                graph_up + 5,
                                float(f(round(graph_low)))
                                - (
                                    float(f(round(graph_low)))
                                    + float(f(round(graph_up)))
                                )
                                / 2
                                - 1,
                                float(f(round(graph_up)))
                                + (
                                    float(f(round(graph_low)))
                                    + float(f(round(graph_up)))
                                )
                                / 2
                                + 1,
                            ]
                        )

                    elif (float(f(graph_low)) == 0) or (float(f(graph_up)) == 0):
                        plt.axis(
                            [
                                graph_low - 5,
                                graph_up + 5,
                                -(graph_up - graph_low) / 2,
                                (graph_up + graph_low) / 2,
                            ]
                        )

            except:
                plt.axis([-7.5, 7.5, -7.5, 7.5])

            plt.legend(loc="lower left")
            plt.grid()

            download = input("Download graph? (y/n) ")

            if download == "y":
                plt.savefig(save_path)
                print("\nDownloaded graph.")

            plt.show()

            return "\nExited graph."

        except:
            plt.close()
            warning("Could not graph function.")
            return "\nExited graph."

    else:
        return "\n[bright_yellow]Exiting Definite Integral Screen ... ... ...[/bright_yellow]"


def improp_int(function: str, low: str, up: str):
    r"""Calculates the improper integral of a function.

    Explanation
    ===========

    Takes a function and a lower and upper bound (can be inf)
    and returns the improper integral from the lower bound to
    the upper bound of the function. Uses Cauchy principal value
    for integrals with infinite bounds and standard integration
    for integrals with singularities.
    """

    function = replace_expr(function)

    check_num(low)

    calc_low = eval(replace_expr(low))

    check_num(up)

    calc_up = eval(replace_expr(up))

    try:
        improper_area = Integral(
            function, (x, calc_low, calc_up)
        ).principal_value()  # Cauchy Principal Value

        if "Integral" in str(improper_area):
            warning("Cannot compute improper integral.\n")
            return ""

        print(
            f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
        )
        print_expr(improper_area)
        improper_area = replace_expr(str(improper_area))

        if history[len(history) - 1] != improper_area:
            history.append(improper_area)

    except ValueError:
        improper_area = integrate(function, (x, calc_low, calc_up))

        if "Integral" in str(improper_area):
            warning("Cannot compute improper integral.\n")
            return ""

        print(
            f"\nCalculated improper integral of [bright_magenta]{function}[/bright_magenta] from {low} to {up}. Final area is:\n"
        )
        print_expr(improper_area)
        improper_area = replace_expr(str(improper_area))

        if history[len(history) - 1] != improper_area:
            history.append(improper_area)

    print()


def double_int(function: str, out_low: str, out_up: str, in_low: str, in_up: str):
    r"""Calculates the double integral of a function over region R.

    Explanation
    ===========

    Takes a function and outer lower and upper and
    inner lower and upper bounds (that define region R)
    and returns the integral of the function over R.
    """

    function = replace_expr(function)

    check_num(out_low)

    out_low = eval(replace_expr(out_low))

    check_num(out_up)

    out_up = eval(replace_expr(out_up))

    in_low = eval(replace_expr(in_low))

    in_up = eval(replace_expr(in_up))

    out_up = eval(replace_expr(str(out_up)))
    out_low = eval(replace_expr(str(out_low)))
    in_up = eval(replace_expr(str(in_up)))
    in_low = eval(replace_expr(str(in_low)))

    result = integrate(function, (y, in_low, in_up), (x, out_low, out_up))

    if "Integral" in str(result):
        warning("Cannot compute integral.\n")
        return ""

    print(
        f"\nDouble integral of [bright_magenta]{function}[/bright_magenta] with inner bounds [bright_cyan]{in_low}[/bright_cyan] and [bright_cyan]{in_up}[/bright_cyan] and outer bounds {out_low} and {out_up} is:\n"
    )
    print_expr(result)
    result = replace_expr(str(result))

    if history[len(history) - 1] != result:
        history.append(result)
    print()


def lim(function: str, value: str):
    r"""Calculates limit of a function at a value.

    Explanation
    ===========

    Takes a function and a value and returns the
    limit of the function at the point.
    """

    function = replace_expr(function)

    check_num(value)

    value = float(eval(replace_expr(value)))

    l = limit(function, x, value)

    if "Limit" in str(l):
        warning("Cannot compute limit.")
        return ""

    if limit(function, x, value, "+") != limit(function, x, value, "-"):
        warning(
            "\nThe limit does not exist (the limit approaching from the right does not equal the limit approaching from the left)."
        )
        print(
            "[bright_red]Use one-side limit calculation to calculate the limit from one direction.[/bright_red]\n"
        )
        return ""

    print(
        f"\nLimit of [bright_magenta]{function}[/bright_magenta] as x approaches {value} is:\n"
    )
    print_expr(l)
    l = replace_expr(str(l))

    if history[len(history) - 1] != l:
        history.append(l)


def side_lim(function: str, value: str, direction: str):
    r"""Calculates limit of a function at a value.

    Explanation
    ===========

    Takes a function and a value and returns the
    limit of the function at the point.
    """
    function = replace_expr(function)

    check_num(value)

    value = float(eval(replace_expr(value)))

    if direction == "left" or direction == "Left":
        direction_sign = "-"

    elif direction == "right" or direction == "Right":
        direction_sign = "+"

    else:
        print()
        error("\nDirection is neither right or left.")
        return ""

    l = limit(function, x, value, dir=direction_sign)

    if "Limit" in str(l):
        warning("\nCannot compute limit.")
        return ""

    print(
        f"\nLimit of [bright_magenta]{function}[/bright_magenta] as x approaches {value} from the [bright_cyan]{direction}[/bright_cyan] is:\n"
    )
    print_expr(l)
    l = replace_expr(str(l))

    if history[len(history) - 1] != l:
        history.append(l)


def eq_solve(mode: str):
    r"""Solves an/a set of equation(s).

    Explanation
    ===========

    Takes a mode to determine the number of equations
    to solve (1 to 3) and returns the values of the variables.
    """

    eq_list = []

    if mode == "1":
        eq1 = input("\nEnter equation: ")
        left = eval(eq1[: eq1.find("=")])
        right = eval(eq1[eq1.find("=") + 1 :])
        eq_set = solve(left - right)
        if len(eq_set) == 0:
            print()
            error("Cannot solve equation")
        else:
            print("\nx:\n")
            for i in range(0, len(eq_set)):
                pprint(eq_set[i])
                sol = replace_expr(str(eq_set[i]))

                if history[len(history) - 1] != sol:
                    history.append(sol)

                print()

    elif mode == "2":
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

    elif mode == "3":
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


################
# Helper funcs #
################


# These assistant functions sanatize inputs and prevent errors


def check_simp(expr: str) -> bool:
    r"""Check for simplification of an expression.

    Explanation
    ===========

    Takes an expression and simplifies/rewrites it
    if possible. Otherwise returns boolean value False.
    """

    try:
        if str(simplify(expr, evaluate=False)) != expr:
            print("\nSimplify/rewrite:\n")
            print_expr(simplify(expr, evaluate=False))

        else:
            return False

    except:
        pass


def check_order(order: str):
    r"""Check if the order of differentiation is valid.

    Explanation
    ===========

    Takes the order of differentiation and gives an error
    if it is invalid.
    """

    if ("." in order) or (order.isnumeric() == False):
        print()
        error("Invalid order of differentiation.")


def check_num(num: str):
    r"""Checks if a string is numerical (valid).

    Explanation
    ===========

    Takes a string and checks if it is numerical valid
    or containing a symbol. It gives an error if it
    is invalid.
    """

    if (
        (num.isnumeric() is False)
        and ("pi" not in num)
        and ("e" not in num)
        and ("E" not in num)
        and ("-" not in num)
        and ("." not in num)
        and ("sqrt" not in num)
        and ("oo" not in num)
        and ("/" not in num)
        and ("x" not in num)
        and ("y" not in num)
        and ("z" not in num)
        and ("t" not in num)
        and ("m" not in num)
    ):
        print()
        error("Invalid integration bound/numerical expression.")


def replace_expr(expr: str) -> str:
    r"""Sanatizes an expression from input.

    Explanation
    ===========

    Takes a string and sanatizes input for calculation.
    """

    expr = expr.replace(" ", "")

    for r in expr_replace:
        expr = expr.replace(*r)

    return expr


def replace_graph(function: str) -> str:
    r"""Replaces an expression (function) for graphing.

    Explanation
    ===========

    Takes a string and sanatizes it for graphing.
    This replacement includes changing prefixes of trig
    functions (since Sympy and Numpy trig functions uses
    different prefixes) and an implementation of csc, sec,
    cot, etc. (since Numpy does not provide them).
    """

    # Graph constant functions

    if "x" not in function:
        function = "0*x+" + function

    for r in graph_replace:
        function = function.replace(*r)

    return function


def replace_calc(expr: str):
    r"""Replaces an expression for calculation.

    Explanation
    ===========

    Takes a string and sanatizes it for simple calculation
    (i.e., calculation without calculus operations). This
    replacement algorithm mainly deals with results that
    should be 0 but instead outputs an answer near 0 (and is
    simply changed to just 0), and complex numbers with an
    imaginary part close to 0 (changed to 0). e.g., exp(i*pi)
    gives a -1.0 * (a really small number). Using replace_calc
    changes it to just -1.0.
    """

    try:
        expr = simplify(expr.replace(" ", "").replace("ss.gamma", "gamma")).evalf()

        if "I" not in str(expr):
            if eval(str(expr)) < 0.000001 and eval(str(expr)) > -0.000001:
                expr = 0

        else:
            expr = str(expr)
            if "+" in str(expr):
                imag = simplify(expr[expr.find("+") + 1 :].rstrip("*I")).evalf()
                if imag < 0.000001 and imag > -0.000001:
                    expr = simplify(expr[: expr.find("+")].rstrip("*I")).evalf()

                else:
                    expr = simplify(expr).evalf()

            else:
                if (
                    simplify(expr.rstrip("*I")).evalf() < 0.000001
                    and simplify(expr.rstrip("*I")).evalf() > -0.000001
                ):
                    expr = 0

                else:
                    expr = simplify(expr).evalf()

        return expr

    except:
        pass


def print_expr(text: str):
    r"""Selects printing method.

    Explanation
    ===========

    Linked to the printing settings option in settings.
    Chooses either normal print or Sympy pretty print.
    Default setting in Sympy pprint.
    """

    printing_methods = {"p": lambda t: pprint(text), "n": lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods["p"](text)


def factorial(x):
    r"""
    Simple implementation of a factorial function,
    using the gamma function as the expansion of
    the factorial.
    """

    return ss.gamma(x + 1)


def nprint(text: str):
    print(text)
    sleep(0.04)


####################
# Driving programs #
####################


def main():
    with Progress(
        SpinnerColumn(finished_text="[bright_green]√[/bright_green]"),
        TextColumn("[bright_green]Handling imports ...[/bright_green]"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        transient=False,
    ) as progress:
        task = progress.add_task("", total=50)

        while not progress.finished:
            progress.update(task, advance=1)
            sleep(randint(2, 5) / 100)

    global x, y, z, t, m, history, expr_replace, graph_replace, save_path

    while True:
        instruct_path = (
            dirname(abspath(__file__))
            + "/texts/main_screen.txt"
            # TODO: make the PATH compatible with Windows
        )
        main = open(instruct_path, mode="r")

        for line in main.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        now = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M:%S")

        print(f"\n(Time now is: {now})\n")
        print("[bold bright_green](Current Screen: Main Screen)[/bold bright_green]\n")
        print("[purple]Enter Command: [/purple]", end="")
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
            print("\n[bright_yellow]Exiting Calc-ULTRA ... ... ...[/bright_yellow]\n")
            break

        else:
            print()
            warning(f'Invalid command:"{cmd}"\n')


"""
If you find this message, type 'hi' in the general discussions - sudoer-Huatao
"""


def derivacalc():
    while True:
        instruct_path = dirname(abspath(__file__)) + "/texts/derivacalc.txt"
        derivacalc = open(instruct_path, mode="r")

        for line in derivacalc.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        print()

        try:
            print(
                "[bold bright_green](Current Screen: DerivaCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                print(
                    "\n[bold bright_green](Current Screen: Derivative Screen)[/bold bright_green]\n"
                )
                function = prompt("Enter a function: ", key_bindings=key_bindings)
                order = prompt(
                    "Enter order of derivative calculation: ", key_bindings=key_bindings
                )

                derive(function, order)

                input("\nPress any key to continue: ")

            elif cmd == "2":
                print(
                    "\n[bold bright_green](Current Screen: Partial Derivative Screen)[/bold bright_green]\n"
                )
                function = prompt(
                    "Enter a function containing x and y or x and y and z: ",
                    key_bindings=key_bindings,
                )
                var = prompt(
                    "Enter variable to differentiate in respect to: ",
                    key_bindings=key_bindings,
                )

                if var != "x" and var != "y" and var != "z":
                    print()
                    error("Invalid variable to differentite in respect to.")

                else:
                    order = prompt(
                        "Enter the order of partial derivative calculation: ",
                        key_bindings=key_bindings,
                    )

                    partial_derive(function, var, order)

                    input("\nPress any key to continue: ")

            elif cmd == "3":
                print(
                    "\n[bold bright_green](Current Screen: Implicit Derivative Screen)[/bold bright_green]\n"
                )
                circ = prompt(
                    "Enter an equation containing x and y:", key_bindings=key_bindings
                )
                order = prompt(
                    "Enter order of implicit derivative calculation: ",
                    key_bindings=key_bindings,
                )

                implicit_derive(circ, order)

                input("\nPress any key to continue: ")

            elif cmd == "4":
                print("\n[bright_yellow]Exiting DerivaCalc ...[/bright_yellow]")
                break

            else:
                print()
                warning(f'Invalid command:"{cmd}"')
                print()

        except:
            print()
            error("An unknown error occured.")
            print("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def intecalc():
    while True:
        instruct_path = dirname((__file__)) + "/texts/intecalc.txt"
        intecalc = open(instruct_path, mode="r")

        for line in intecalc.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        print()

        try:
            print(
                "[bold bright_green](Current Screen: InteCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                print(
                    "\n[bold bright_green](Current Screen: Antiderivative Screen)[/bold bright_green]\n"
                )
                function = prompt("Enter a function: ", key_bindings=key_bindings)

                antiderive(function)

                input("\nPress any key to continue: ")

            elif cmd == "2":
                print(
                    "\n[bold bright_green](Current Screen: Definite Integral Screen)[/bold bright_green]\n"
                )
                function = prompt("Enter a function: ", key_bindings=key_bindings)
                lower_bound = prompt(
                    "\nEnter the lower bound: ", key_bindings=key_bindings
                )
                upper_bound = prompt(
                    "Enter the upper bound: ", key_bindings=key_bindings
                )

                print(def_int(function, lower_bound, upper_bound))

                input("\nPress any key to continue: ")

            elif cmd == "3":
                print(
                    "\n[bold bright_green](Current Screen: Improper Integral Screen)[/bold bright_green]\n"
                )
                function = prompt("Enter a function: ", key_bindings=key_bindings)
                lower_bound = prompt(
                    "\nEnter the lower bound: ", key_bindings=key_bindings
                )
                upper_bound = prompt(
                    "Enter the upper bound: ", key_bindings=key_bindings
                )

                improp_int(function, lower_bound, upper_bound)

                input("\nPress any key to continue: ")

            elif cmd == "4":
                print(
                    "\n[bold bright_green](Current Screen: Double Integral Screen)[/bold bright_green]\n"
                )
                function = prompt("Enter a function: ", key_bindings=key_bindings)
                outer_low = prompt(
                    "\nEnter the lower outer bound: ", key_bindings=key_bindings
                )
                outer_up = prompt(
                    "Enter the upper outer bound: ", key_bindings=key_bindings
                )
                inner_low = prompt(
                    "\nEnter the lower inner bound: ", key_bindings=key_bindings
                )
                inner_up = prompt(
                    "Enter the upper inner bound: ", key_bindings=key_bindings
                )

                double_int(function, outer_low, outer_up, inner_low, inner_up)

                input("\nPress any key to continue: ")

            elif cmd == "5":
                print("\n[bright_yellow]Exiting InteCalc ...[/bright_yellow]")
                break

            else:
                print()
                warning(f'Invalid command: "{cmd}"')

        except:
            print()
            error("An unknown error occured.")
            print("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def limcalc():
    while True:
        instruct_path = dirname(abspath(__file__)) + "/texts/limcalc.txt"
        limcalc = open(instruct_path, mode="r")

        for line in limcalc.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        print()

        try:
            print(
                "\n[bold bright_green](Current Screen: LimCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                print(
                    "\n[bold bright_green](Current screen: Limit Screen)[/bold bright_green]\n"
                )
                expr = prompt("Enter an expression: ", key_bindings=key_bindings)
                value = prompt("Enter point of evaluation: ", key_bindings=key_bindings)

                lim(expr, value)

                input("\nPress any key to continue: ")

            elif cmd == "2":
                print(
                    "\n[bold bright_green](Current screen: One-sided Limit Screen)[/bold bright_green]\n"
                )
                expr = prompt("Enter an expression: ", key_bindings=key_bindings)
                value = prompt("Enter point of evaluation: ", key_bindings=key_bindings)
                direction = input("Enter direction of limit ('left' or 'right'): ")

                side_lim(expr, value, direction)

                input("\nPress any key to continue: ")

            elif cmd == "3":
                print("\n[bright_yellow]Exiting LimCalc ...[/bright_yellow]")
                break

            else:
                print()
                warning(f'Invalid command: "{cmd}"')

        except:
            print()
            error("An unknown error occured.")
            print("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def algcalc():
    while True:
        instruct_path = dirname(abspath(__file__)) + "/texts/algcalc.txt"
        algcalc = open(instruct_path, mode="r")

        for line in algcalc.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        print()
        try:
            print(
                "\n[bold bright_green](Current Screen: AlgCalc Main Screen)[/bold bright_green]\n"
            )
            print("[purple]Enter Command: [/purple]", end="")
            cmd = input()

            if cmd == "1":
                print(
                    "\n[bold bright_green](Current screen: Equation Solver Screen)[/bold bright_green]\n"
                )
                mode = ""

                while mode != "q":
                    print(
                        'Enter mode: 1 for one set equation, 2 for two, and 3 for three ("q" to quit): ',
                        end="",
                    )
                    mode = input()
                    eq_solve(mode)

            elif cmd == "2":
                print(
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
                    expr = prompt(key_bindings=key_bindings)
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
                            ).replace(" ", "")
                            result = eval(replace_graph(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]", end="")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        error(f'Could not parse: "{expr}"\n')

                print()

            elif cmd == "3":
                print(
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
                    expr = prompt(key_bindings=key_bindings)
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
                            ).replace(" ", "")
                            result = eval(replace_graph(expr))
                            print("\n[bright_yellow]Result: [/bright_yellow]\n")
                            pprint(result)
                            print()

                        else:
                            break
                    except:
                        print()
                        error(f'Could not parse: "{expr}"\n')

            elif cmd == "4":
                print("\n[bright_yellow]Exiting AlgCalc ...[/bright_yellow]")
                break

            else:
                print()
                warning(f'Invalid command: "{cmd}"')

        except:
            print()
            error("An unknown error occured.")
            print("[bold bright_red]Check if your input is valid.[/bold bright_red]\n")


def settings():
    while True:
        settings_path = dirname(abspath(__file__)) + "/texts/settings.txt"
        settings = open(settings_path, mode="r")

        for line in settings.readlines():
            line = line.rstrip()
            if line.startswith(" "):
                nprint(f"[gold1]{line}[/gold1]")
            else:
                nprint(line)

        print("\n[bold green](Current Screen: Settings Screen)[/bold green]\n")
        print("[purple]Enter Command: [/purple]", end="")
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

        elif cmd == "style":
            print(
                "\n[bold bright_green](Current Screen: Graph Style Settings Screen)[/bold bright_green]\n"
            )

            print("You currently have these available styles:\n")

            style_list = plt.style.available
            styles = "  ,  ".join(style_list)
            nprint(f"[bright_yellow]{styles}[/bright_yellow]")

            style = input('\nChoose a style to apply (type "default" to reset): ')

            plt.style.use(style)

            print(f'Graph style set to "{style}".')

        elif cmd == "path":
            print(
                "\n[bold bright_green](Current Screen: Graph Save Path Screen)[/bold bright_green]\n"
            )

            global save_path

            print(f"The current path where graphs will be saved to is: {save_path}.\n")

            path = input(
                'Enter a new absolute path for the graph to be saved in (e.g. /Users/spam/eggs), or exit ("q"): '
            )

            if path != "q":
                if exists(path):
                    save_path = path + "/fig.pdf"
                    print(f"\nPath saved to: {save_path}.")

                else:
                    warning("\nThe path you specificied does not exist.")

        elif cmd == "exit":
            print("\n[bright_yellow]Exiting settings ... ... ...[/bright_yellow]")
            break

        else:
            print()
            warning(f'Invalid command:"{cmd}"')


if __name__ != "__main__":
    main()
# You've reached the end of the file!
