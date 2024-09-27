===========================
Tutorial for sensitivity analysis of the Vanderpol equation
===========================

Introduction
------------
This tutorial demonstrates how to perform sensitivity analysis of the Vanderpol equation using the library provided by this repository.

The Vanderpol equation is a second-order nonlinear differential equation, and is given by:

$$
\begin{bmatrix}
	\dot{u}_1\\
	\dot{u}_2
\end{bmatrix}
=
\begin{bmatrix}
	u_2\\
	\mu\left((1-u_1^2)u_2-u_1\right)
\end{bmatrix},
$$

Installation
------------

To install the library, run:

.. code-block:: bash

   git clone https://github.com/username/repo.git
   cd repo
   make install

Usage
-----

To use the library, you can do the following:

.. code-block:: cpp

   #include "my_library.h"

   int main() {
       // Your code here
   }
