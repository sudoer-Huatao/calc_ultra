# calc-ultra

[![GPLv3 License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/license/mit/) [![Version](https://img.shields.io/badge/Version-1.3.3-blue.svg)](https://github.com/sudoer-Huatao/Calc-ULTRA_Calculus-Calculator)

> **Calculus made easy**

(Turn on dark mode for a better aesthetic)

The Calc-ULTRA calculus calculator, but as a module!

- Little Python background knowledge needed!

Supports:

- Derivatives
- Partials
- Implicit differentiation
- Antiderivatives
- Definite integrals
- Improper integrals
- Double integrals
- Solving (sets) of equation(s)
- Vector/matrix operations
- **A perfect interface to do calculations!**  

## Note

This is the module package of the Calc-ULTRA calculator. For the Python script of this package, visit <https://github.com/sudoer-Huatao/Calc-ULTRA> (**unmaintained**).

## Installation and Running

> Run the calculus calculator with a single line of code

Command line: `pip3 install calc-ultra`.
Due to Python import identifiers restrictions, please import Calc-ULTRA as "calc_ultra" and not "calc-ultra" when you need to use the calculator.

To run the calculator, import Calc-ULTRA as `calc_ultra` like so:

`from calc_ultra import main`

Make sure you have the latest version installed. To update calc-ultra, run `pip3 install --upgrade calc-ultra`.

## Requirements

This program requires the `sympy`,  `numpy`, `rich`, `matplotlib`, and `scipy` modules installed. Other required modules are built in to most Python IDEs.

## Warnings

### Function limitations

Due to limitations of the SymPy module, **some functions cannot be integrated**. The Error Function `erf(x)` can be integrated in both indefinite integral and definite integral calculation, but the Absolute Value and Factorial functions are only available to definite integral calculations. Integration of composed functions are also limited due to SymPy limitations. While some composed functions work, others don't. 😟

## Test PYPI

Previous test versions of this project are on Test PYPI. View on <https://test.pypi.org/project/calc-ultra/>.

## Acknowledgements

> Without them, this would be impossible

A general thank-you to all GitHub users who gave feedback and/or starred this repository. ⭐️
And... a SPECIAL THANK-YOU to @Haobot for troubleshooting and feedback! 👍❤️

This program was made using SymPy and Scipy for calculation and Matplotlib and NumPy for graphing.

## Gallery

DerivaCalc derivative with graph demo:
![derivacalc_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/22f7f674-81b4-4e5f-9bff-586d94c1976f "derivacalc_demo")

InteCalc antiderivative with graph demo:
![intecalc_demo_1](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/1035a5af-cf04-46c1-b2e3-ec5eff08c58d "intecalc_demo_1")

InteCalc definite integral with graph demo:
![intecalc_demo_2](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/42a127af-17f6-4608-a6aa-9029cf00e973 "intecalc_demo_2")

LimCalc limit demo:
![limcalc_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/bb08db05-9b1b-4204-9d98-d4a3baa2167d "limcalc_demo")

AlgCalc equation solver demo:
![algcalc_demo_1](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/7d2bdad0-572c-4eec-8644-2837a0689154 "algcalc_demo_1")

AlgCalc vector operation demo:
![algcalc_demo_2](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/28c5143b-93fc-4165-b8d9-572d1ab00da8 "algcalc_demo_2")

AlgCalc matrix operation demo:
![algcalc_demo_3](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/b6d11f38-b5b8-42ed-b4b3-664c65207664 "algcalc_demo_3")

## License

This project is licensed under the terms of the MIT license.
