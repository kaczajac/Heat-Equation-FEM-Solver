# Heat Equation Solver

A python program that solves the heat equation using finite element method.

## Problem Statement

Find the function _**u**_ : [0,2] -> R satisfying the equation:
  
$$
\begin{equation}
-\frac{d}{dx} \left( k(x) \frac{du(x)}{dx} \right) = 100x \quad x \in [0, 2] \\
\\
\end{equation}
$$

With boundary conditions:

$$
\begin{equation}
u(2) = 0 \\
\\
\end{equation}
$$

$$
\begin{equation}
\frac{du(0)}{dx} + u(0) = 20 \\
\\
\end{equation}
$$

Where function _**k**_ is defined as follows:

$$
\begin{equation}
k(x) =
\begin{cases}
1 & \text{for } x \in [0, 1] \\
2x & \text{for } x \in (1, 2]
\end{cases}
\\  
\end{equation}
$$

## Contents
- _**main.py**_ -> main module
- _**helpers.py**_ -> contains helper functions for _main.py_ module 
- _**weak-formulation.pdf**_ -> contains a step-by-step procedure for finding the weak formulation of the given equation.

## How to run
Make sure _helpers.py_ module is in the same directory as _main.py_. You will also need to install _matplotlib_ library - it will be needed
for plotting the final result.

After that you can simply execute _main.py_ module:
```
$ python main.py
```
The program then prompts the user to specify the number of elements for FEM procedure.
