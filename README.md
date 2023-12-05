### (Note: This repo is the calc_ultra module. For the Python script, visit https://github.com/sudoer-Huatao/Calc-ULTRA)

# Calc-ULTRA: Calculus Calculator

[![GPLv3 License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/license/mit/) [![Version](https://img.shields.io/badge/Version-1.0.0-blue.svg)](https://github.com/sudoer-Huatao/Calc-ULTRA_Calculus-Calculator)

Calc-ULTRA, but a module!

- No Python background knowledge needed.
- Supports derivatives, antiderivatives, and definite integrals - with a graph!
- (Please FULL SCREEN your console/terminal for better experience!)
<img width="1710" alt="截屏2023-11-28 15 08 37" src="https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/bf1e55dd-cb70-46cf-bdc9-0bf58a7886c6">

## How to run

To run it, simply use pip to install. Then, import main from the calc_ultra module. That should run the calculator. Demo:

https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/bfa34921-9d15-4160-af36-2cb9ccd912b9


## Requirements

This program requires `sympy`,  `numpy`, `matplotlib`, `datetime`, `math`, `logging`, `warnings` and `os` modules. `datetime`, `math`, `logging`, `warnings` and `os` are built-in to most Python editors. The rest can be installed with the command `pip install MODULE_NAME`, though the required modules should be automatically downloaded.

<img width="1710" alt="截屏2023-11-28 15 11 12" src="https://github.com/sudoer-Huatao/calc_ultra/assets/135504586/20e1b072-d2dd-4610-9378-95fd74b80a84">

# Warnings!!!

## Function limitations:

Due to limitations of the SymPy module, SOME FUNCTIONS CAN NOT BE INTEGRATED. The Error Function erf(x) can be integrated in both indefinite integral and definite integral calculation, but the Absolute Value and Factorial functions are only available to definite integral calculations. Also, the factorial function cannot be graphed properly. Integration of composed functions are also limited due to SymPy limitations. While some composed functions work, others don't. 😟

## Test PYPI

The test version of this project is on Test PYPI. View on https://test.pypi.org/project/calc-ultra/.

# Acknowledgements

A general thank-you to all GitHub users who gave feedback and/or starred this repository. ⭐️
And... a SPECIAL THANK-YOU to @Haobot for troubleshooting and feedback! 👍❤️

This program was made using SymPy and SciPy for calculation and Matplotlib and NumPy for graphing.

# License

This project is licensed under the terms of the MIT license.
