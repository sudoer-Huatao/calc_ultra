# calc-ultra

[![GPLv3 License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/license/mit/) [![Version](https://img.shields.io/badge/Version-1.3.1-blue.svg)](https://github.com/sudoer-Huatao/Calc-ULTRA_Calculus-Calculator)

> **Calculus made easy**

The Calc-ULTRA calculus calculator, but as a module!

- Little Python background knowledge needed!

Supports:

- Derivatives
- Partials
- Implicit differentiation
- Antiderivatives
- Definite integrals
- Improper integrals

**NEW SINCE V1.2.1:** Graphs for differentiation and integrals are supported!

## Note

This is the module package of the Calc-ULTRA calculator. For the Python script of this package, visit <https://github.com/sudoer-Huatao/Calc-ULTRA>

## Installation and Running

> Run the calculus calculator with a single line of code

Command line: `pip3 install calc-ultra`.
Due to Python import identifiers restrictions, please import Calc-ULTRA as "calc_ultra" and not "calc-ultra" when you need to use the calculator.

To run the calculator, import Calc-ULTRA as `calc_ultra` like so:

`from calc_ultra import main`

Make sure you have the latest version installed. To update calc-ultra, run `pip3 install --upgrade calc-ultra`.

## Requirements

This program requires the `sympy`,  `numpy`, `rich`, and `matplotlib` modules installed. Other required modules are built in to most Python IDEs.

## Warnings

### Function limitations

Due to limitations of the SymPy module, **some functions cannot be integrated**. The Error Function `erf(x)` can be integrated in both indefinite integral and definite integral calculation, but the Absolute Value and Factorial functions are only available to definite integral calculations. Also, the factorial function cannot be graphed properly. Integration of composed functions are also limited due to SymPy limitations. While some composed functions work, others don't. 😟

## Test PYPI

Previous test versions of this project are on Test PYPI. View on <https://test.pypi.org/project/calc-ultra/>.

## Acknowledgements

> Without them, this would be impossible

A general thank-you to all GitHub users who gave feedback and/or starred this repository. ⭐️
And... a SPECIAL THANK-YOU to @Haobot for troubleshooting and feedback! 👍❤️

This program was made using SymPy for calculation and Matplotlib and NumPy for graphing.

## Gallery

DerivaCalc derivative with graph demo:
![derivacalc_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/2cac3688-d1da-4b73-b121-b6a2c7cf471c "derivacalc_demo")

InteCalc antiderivative with graph demo:
![intecalc_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/bcfea09c-078f-4344-9ff0-91a414614244 "intecalc_demo")

InteCalc definite integral with graph demo:
![intecalc_graph_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/573e1b76-98f8-492a-845a-721f8d7ccd63 "intecalc_graph_demo")

LimCalc limit demo:
![limcalc_demo](https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/c1069e1c-27ef-4b7a-b68b-49432ed31e5f "limcalc_demo")

## License

This project is licensed under the terms of the MIT license.
