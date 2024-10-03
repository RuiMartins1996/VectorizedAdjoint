Tutorial for sensitivity analysis of the van der Pol oscillator
==================================================================

Introduction
------------
The Van Der Pol oscillator is a second-order nonlinear ordinary differential equation,

.. math::
    \ddot{x} - \mu (1 - x^2) \dot{x} + \mu x = 0


In Wikipedia, the 0th order term is not multiplied by :math:`\mu`, but in this tutorial we will consider the same equation used in the PETSc TSAdjoint example 
`ex20adj.c <https://petsc.org/release/src/ts/tutorials/ex20adj.c.html>`_
(notice how RHSFunction() is defined).

This can be rewritten as a system of first-order ordinary differential equations,

.. math::
    \begin{bmatrix}
	    \dot{x}\\
	    \dot{v}
    \end{bmatrix}
    =
    \begin{bmatrix}
	    v\\
        -\mu(x -(1-x^2) v)
    \end{bmatrix}.

where :math:`\mu` is a damping coefficient.

Using an adaptive Stepper
--------------------------------
In this example, we will use an adaptive stepper, the Runge-Kutta-Fehlberg 7(8) error stepper available in :code:`boost::numeric::odeint`. We additionally define a constant :code:`tol` to be used as the absolute and relative tolerances of the stepper.

.. code-block:: cpp
    
    const double tol = 1e-5;
    
    typedef odeint::runge_kutta_fehlberg78<std::vector<double>> stepper_type;

To perform the forward pass, we must pass an instance of a controlled stepper to the :code:`runge_kutta()` function. The controlled stepper is created by passing the stepper type and the error tolerances to the :code:`odeint::make_controlled` function.

.. code-block:: cpp

    // Absolute and relative tolerance
    double AbsTol = tol, RelTol = tol;

    auto ctrlStepper = odeint::make_controlled<stepper_type>(AbsTol, RelTol);

The way adaptivity works and more details on the stepperscan be found `here <https://live.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html>`_.




Performing sensitivity analysis
--------------------------------
We now wish to find the sensitivities of the solution at time :math:`t = t_f` with respect to the initial conditions :math:`\left[x(t_i),v(t_i)\right]^T` and the damping coefficient :math:`\mu`. The objective functions are the final values of the state vector 

.. math::
    \begin{bmatrix}
        x(t_f)\\
        v(t_f)
    \end{bmatrix}

We create an instance of the :code:`Driver` class,

.. code-block:: cpp
   
   // System size
   int Nin = 2;
   // Number of cost functions 
   int Nout = 2;
   // Number of parameters
   int Npar = 1;

   Driver driver(Nin, Nout, Npar);

Forward Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will consider the following initial conditions and damping coefficient,

.. code-block:: cpp

    // Time interval and time step
    double ti = 0.0 ,tf = 5e-1,dt = 0.001;

    std::vector<double> mu = {1e3};
    std::vector<double> x0 = {2.0, -2.0 / 3.0 + 10.0 / (81.0 * mu[0]) - 292.0 / (2187.0 * mu[0] * mu[0])};

We then solve the ODE by doing a forward pass,

.. code-block:: 

    // Create the system function instance
    auto vdp = VanDerPol();

    // Forward pass
    size_t numSteps = runge_kutta(ctrlStepper, vdp, x0, mu, ti, tf, dt, driver);

Reverse Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We set the partial derivatives of the cost function with respect to the state vector solution of the ODE (:math:`\left[x(t_f),v(t_f)\right]^T`) and with respect to the parameter :math:`\mu`,

.. code-block:: cpp

    auto lambda = std::vector(Nout, std::vector<double>(Nin));
    lambda[0][0] = 1.0 * r0[0]; // dE(t_f)dx(t_f) = k * x(t_f)
    lambda[0][1] = 1.0 * r0[1]; // dE(t_f)dv(t_f) = m * v(t_f)

    auto muadj = std::vector(Nout, std::vector<double>(Npar));
    muadj[0][0] = 0.0; // dE(t_f)d\mu = 0.0

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, muadj);

Additionally, we need to inform the driver of the chosen stepper'same Butcher Tableau,

.. code-block:: 

    constructDriverButcherTableau(driver, stepper);

and to record the rhs function with automatic differentiation,

.. code-block:: 

    recordDriverRHSFunction(driver, hm);


we perform a reverse pass by doing:

.. code-block::

    // Reverse pass to obtain the adjoints of the cost functions
    backpropagation::adjointSolve(driver, mu);

The sensitivities are stored in :code:`lambda` and :code:`muadj`.

TODO