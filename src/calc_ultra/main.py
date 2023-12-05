from sympy import *
from sympy.core.numbers import pi, E
from math import floor, ceil
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.patches import Polygon
import logging
import os
import warnings


# This program requires the Sympy, NumPy, and MatPlotLib modules installed!

def derive(function, order):
    if '^' in function:
        function = function.replace('^', '**')

    if check_order(order) is True:
        df = diff(function, x, order)

        print(f'\nDerivative of {function} with order {order} is:\n')
        print_expr(df)

        if check_simp(df) is True:
            print('\nSimplify/rewrite:\n')
            print_expr(simplify(df, evaluate=False))


def partial_derive(function, var, order):
    if '^' in function:
        function = function.replace('^', '**')

    if check_order(order) is True:
        df = diff(function, var, order)

        print(f'\nPartial derivative of {function} in respect to {var} of order {order} is:\n')
        print_expr(df)

        if check_simp(df) is True:
            print('\nSimplify/rewrite:\n')
            print_expr(simplify(df, evaluate=False))


def implicit_derive(circ, order):
    if '^' in circ:
        circ = circ.replace('^', '**')

    if check_order(order) is True:
        df = idiff(eval(circ), y, x, int(order))

        print(f'\nDerivative of {circ} with order {order} is:\n')
        print_expr(df)

        if str(simplify(df, evaluate=False)) != str(df):
            print('\nSimplify/rewrite:\n')
            print_expr(simplify(df, evaluate=False))


def antiderive(function):
    if 'pi' in function:
        function = function.replace('pi', str(pi))

    if '^' in function:
        function = function.replace('^', '**')

    F = Integral(function, x).doit()

    if 'Integral' in str(F):
        logging.warn('Cannot compute integral.\n')

    else:
        print('\nAntiderivative is:\n')
        print_expr(F)

        if check_simp(F) is True:
            print('\nSimplify/rewrite:\n')
            print_expr(simplify(F, evaluate=False))

        print("\nDon't forget to add a constant!\n")


def definite_integrate(function, low, up):
    x = Symbol('x')

    if '^' in function:
        function = function.replace('^', '**')

    if check_bound(low) is False:
        return ''

    low = float(eval(replace_bound(low)))

    if check_bound(up) is False:
        return ''

    up = float(eval(replace_bound(up)))

    result = integrate(function, (x, low, up)).evalf()

    if (
        ('1/x' in function or function == 'x^-1')
        and (low <= 0 or low <= low + 1)
        or (str(result) == 'nan')
        or ('I' in str(result))
    ):
        logging.warn('Cannot compute integral because integral does not converge.')

    else:
        print(f'\nCalculated integral of {function} from {low} to {up}. Final area is:\n')
        print_expr(result)
        print('\nShow graph of area? (y/n)')
        show = input('(Exit the graph window when you are finished to continue) ')

        if show == 'y':
            try:
                print('\nLoading graph. Might take some time on first startup ...')

                x = np.linspace((-up - 8), (up + 8), 200000)

                if 'ln' in function or 'log' in function:
                    x = np.linspace(int(floor(low)) + 1, int(ceil(up)) + 8, 200000,)

                y = [g(a) for a in x]
                fig, ax = plt.subplots()

                title = 'Shaded area beneath function'
                plt.title(title)
                plt.xlabel('x', weight='bold')
                plt.ylabel('y', rotation=0, weight='bold')
                plt.plot(x, y, color='red')

                try:
                    if graph_option == 'f':
                        plt.axis([-7.5, 7.5, -7.5, 7.5])

                    elif graph_option == 'a':
                        if (float(g(low)) != 0) and (float(g(up)) != 0):
                            plt.axis(
                                [
                                    low - 5,
                                    up + 5,
                                    float(g(round(low)))
                                    - (
                                        float(g(round(low)))
                                        + float(g(round(up)))
                                    )
                                    / 2
                                    - 1,
                                    float(g(round(up)))
                                    + (
                                        float(g(round(low)))
                                        + float(g(round(up)))
                                    )
                                    / 2
                                    + 1,
                                ]
                            )
                        elif (float(g(low)) == 0) or (float(g(up)) == 0):
                            plt.axis(
                                [
                                    low - 5,
                                    up + 5,
                                    -(up - low) / 2,
                                    (up + low) / 2,
                                ]
                            )

                except:
                    plt.axis([-7.5, 7.5, -7.5, 7.5])

                plt.grid()

                ix = np.linspace(low, up)
                iy = [g(i) for i in ix]
                verts = [(low, 0)] + list(zip(ix, iy)) + [(up, 0)]
                poly = Polygon(verts, facecolor='blue')
                ax.add_patch(poly)

                plt.show()
                return '\nExited graph.'

            except:
                logging.warn('Could not graph function.')

        else:
            return '\nExiting Definite Integral Screen ... ... ...\n'


def improper_integrate(function, low, up):
    if 'pi' in function:
        function = function.replace('pi', str(pi))

    if '^' in function:
        function = function.replace('^', '**')
    else:
        str(function)

    if check_bound(low) is False:
        return ''

    if 'oo' in low:
        low = eval(low)
    else:
        low = float(eval(replace_bound(low)))

    if check_bound(up) is False:
        return ''

    if 'oo' in up:
        up = eval(up)
    else:
        up = float(eval(replace_bound(up)))

    try:
        improper_area = Integral(function, (x, low, up)).principal_value()

        print(f'Calculated improper integral of {function} from {low} to {up}. Final area is:\n')
        print_expr(improper_area)
        print()

    except ValueError:
        logging.warn('ValueError: Singularity while computing improper integral.\n')


def normal_limit(expr, value):
    print('\n(Current screen: Limit Screen)\n')

    if 'pi' in expr:
        expr = expr.replace('pi', str(pi))

    if '^' in expr:
        expr = expr.replace('^', '**')

    if check_bound(value) is False:
        return ''

    value = float(eval(replace_bound(value)))

    l = limit(expr, x, value)

    if 'Limit' in str(l):
        logging.warn('Cannot compute limit.')

    else:
        print(f'\nLimit of {expr} as x approaches {value} is:\n')
        print_expr(l)
        if check_simp(l) is True:
            print('\nSimplify/rewrite:\n')
            print_expr(simplify(l, evaluate=False))


def one_side_limit(expr, value, direction):
    print('\n(Current screen: One-sided Limit Screen)\n')

    if 'pi' in expr:
        expr = expr.replace('pi', str(pi))

    if '^' in expr:
        expr = expr.replace('^', '**')

    if check_bound(value) is False:
        return ''

    value = float(eval(replace_bound(value)))

    if direction == 'left':
        direction_sign = '-'

    elif direction == 'right':
        direction_sign = '+'

    else:
        logging.error('\nTypeError: Direction is neither right or left.')
        return ''

    l = limit(expr, x, value, dir=direction_sign)

    if 'Limit' in str(l):
        logging.warn('\nCannot compute limit.')

    else:
        print(f'\nLimit of {expr} as x approaches {value} from the {direction} is:\n')
        print_expr(l)


def check_simp(expr):
    if str(simplify(expr, evaluate=False)) != str(expr):
        return True
    else:
        return False


def check_order(order):
    if ('.' in order) or (order.isnumeric() == False) or (int(order) <= 0):
        logging.error('OrderError: Order of derivative calculation is not a valid number.')
        return False
    else:
        return True


def check_bound(bound):
    if (
        (bound.isnumeric() is False)
        and ('pi' not in bound)
        and ('e' not in bound)
        and ('-' not in bound)
        and ('.' not in bound)
        and ('sqrt' not in bound)
        and ('oo' not in bound)
    ):
        logging.error('TypeError: Integration bound is a not a number.')
        return False
    else:
        return True


def replace_bound(bound):
    if 'pi' in bound:
        bound = bound.replace('pi', str(pi))
    if 'e' in bound:
        bound = bound.replace('e', str(E))
    return bound


def g(x):
    final = eval(di_function)
    return final


def settings():
    settings_path = os.path.dirname(os.path.abspath(__file__)) + '/texts/settings.txt'
    settings = open(settings_path, mode='r')

    for line in settings.readlines():
        line = line.rstrip()
        print(line)

    while True:

        print('\n(Current Screen: Settings Screen)\n')
        cmd = input('Enter Command: ')

        if cmd == 'print':
            print('\n(Current Screen: Print Settings Screen)\n')

            global print_option
            print_option = input('Set print mode: "p" (Sympy Pretty Print) or "n" (Normal Print): ')
            print(f'\nPrinting mode set to: "{print_option}"')

        elif cmd == 'graph':
            print('\n(Current Screen: Graph Settings Screen)\n')

            global graph_option
            graph_option = input('Set graph mode: "f" (Fixed graph view) or "a" (Adjusted graph view): ')
            print(f'\nGraph mode set to: "{graph_option}"')

        elif cmd == 'date':
            print('\n(Current Screen: Date Settings Screen)\n')

            global date_option
            date_option = input('Set date mode: "1" (YY/MM/DD) or "2" (YY/MM/DD/HH/MM/SS): ')
            print(f'\nDate mode set to: "{date_option}"')

        elif cmd == 'exit':
            print('\nExiting settings ... ... ...')
            break

        else:
            logging.warn(f'Invalid command:"{cmd}"')


def print_expr(text):
    printing_methods = {'p': lambda t: pprint(text), 'n': lambda t: print(text)}

    try:
        printing_methods[print_option](text)

    except NameError:
        printing_methods['p'](text)


def main():
    global x, y, z
    x, y, z = symbols('x,y,z')

    instruct_path = os.path.dirname(os.path.abspath(__file__)) + '/texts/main_screen.txt'
    main = open(instruct_path, mode='r')
    for line in main.readlines():
        line = line.rstrip()
        print(line)

    warnings.filterwarnings('ignore')

    while True:
        try:
            if date_option == '1':
                now = (datetime.datetime.now()).strftime('%Y/%m/%d')
            elif date_option == '2':
                now = (datetime.datetime.now()).strftime('%Y/%m/%d %H:%M:%S')
        except:
            now = (datetime.datetime.now()).strftime('%Y/%m/%d %H:%M:%S')

        print(f'\n(Time now is: {now})')
        print('(Current Screen: Main Screen)\n')
        cmd = input('Enter Command: ')

        if cmd == '1':
            derivacalc()

        elif cmd == '2':
            intecalc()

        elif cmd == '3':
            limcalc()

        elif cmd == '4':
            settings()

        elif cmd == '5':
            print('\nExiting Calc-ULTRA ... ... ...\n')
            break

        else:
            logging.warn(f'Invalid command:"{cmd}"\n')


'''
If you find this message, type 'hi' in the general discussions - sudoer-Huatao
'''


def derivacalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + '/texts/derivacalc_instructs.txt'
    derivacalc = open(instruct_path, mode='r')
    for line in derivacalc.readlines():
        line = line.rstrip()
        print(line)

    while True:
        print('\n(Current Screen: DerivaCalc Main Screen)\n')
        cmd = input('Enter Command: ')

        if cmd == '1':
            print('\n(Current Screen: Derivative Screen)\n')
            function = input('Enter a function: ')
            order = input('Enter order of derivative calculation: ')
            derive(function, order)

        elif cmd == '2':
            print('\n(Current Screen: Partial Derivative Screen)\n')
            function = input('Enter a function containing x and y or x and y and z: ')
            var = input('Enter variable to differentiate in respect to: ')
            if var != 'x' and var != 'y' and var !='z':
                logging.error('Variable to differentite in respect to is invalid.')
            else:
                order = input('Enter the order of partial derivative calculation: ')
                partial_derive(function, var, order)

        elif cmd == '3':
            print('\n(Current Screen: Implicit Derivative Screen)\n')
            circ = input('Enter the left side of an equation containing x and y: (right side default as 0) ')
            order = input('Enter order of implicit derivative calculation: ')
            implicit_derive(circ, order)

        elif cmd == '4':
            print('\nExiting DerivaCalc ... ... ...')
            break

        else:
            logging.warn(f'Invalid command:"{cmd}"')


def intecalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + '/texts/intecalc_instructs.txt'
    intecalc = open(instruct_path, mode='r')
    for line in intecalc.readlines():
        line = line.rstrip()
        print(line)

    while True:
        print('(Current Screen: InteCalc Main Screen)\n')
        cmd = input('Enter Command: ')

        if cmd == '1':
            print('\n(Current Screen: Antiderivative Screen)\n')
            function = input('Enter a function: ')
            antiderive(function)

        elif cmd == '2':
            print('\n(Current Screen: Definite Integral Screen)\n')
            global di_function
            di_function = input('Enter a function: ')
            lower_bound = input('\nEnter the lower bound: ')
            upper_bound = input('Enter the upper bound: ')
            print(definite_integrate(di_function, lower_bound, upper_bound))

        elif cmd == '3':
            print('\n(Current Screen: Improper Integral Screen)\n')
            function = input('Enter a function: ')
            lower_bound = input('\nEnter the lower bound: ')
            upper_bound = input('Enter the upper bound: ')
            improper_integrate(function, lower_bound, upper_bound)

        elif cmd == '4':
            print('\nExiting InteCalc ... ... ...')
            break

        else:
            logging.warn(f'Invalid command: "{cmd}"')


def limcalc():
    instruct_path = os.path.dirname(os.path.abspath(__file__)) + '/texts/limcalc_instructs.txt'
    limcalc = open(instruct_path, mode='r')
    for line in limcalc.readlines():
        line = line.rstrip()
        print(line)

    while True:
        try:
            print('\n(Current Screen: LimCalc Main Screen)\n')
            cmd = input('Enter Command: ')

            if cmd == '1':
                expr = input('Enter an expression: ')
                value = input('Enter point of evaluation: ')
                normal_limit(expr, value)

            elif cmd == '2':
                expr = input('Enter an expression: ')
                value = input('Enter point of evaluation: ')
                direction = input("Enter direction of limit ('left' or 'right'): ")
                one_side_limit(expr, value, direction)

            elif cmd == '3':
                print('\nExiting LimCalc ... ... ...')
                break

            else:
                logging.warn(f'Invalid command: "{cmd}"')

        except:
            logging.error('UnknownError: An unknown error occured.')


main()
