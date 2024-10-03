Van Der Pol oscillator
==================================================================

Objective
------------
To find how the solution at :math:`t = t_f` varies with respect to the initial conditions and the damping coefficient.

Introduction
------------
The Van Der Pol oscillator is a second-order nonlinear ordinary differential equation:

.. math::
    \ddot{x} - \mu (1 - x^2) \dot{x} + \mu x = 0

In Wikipedia, the 0th order term is not multiplied by :math:`\mu`, but in this tutorial we will consider the same equation used in the PETSc TSAdjoint example 
`ex20adj.c <https://petsc.org/release/src/ts/tutorials/ex20adj.c.html>`_
(notice how RHSFunction() is defined).

This equation can be cast into the following system of ODEs:

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
In this example, we use an adaptive stepper, the Runge-Kutta-Fehlberg 7(8) error stepper, available in :code:`boost::numeric::odeint`.

.. code-block:: cpp

    using namespace boost::numeric::odeint;
    typedef odeint::runge_kutta_fehlberg78<std::vector<double>> stepper_type;

To perform the forward pass, we must pass an instance of a controlled stepper to the :code:`runge_kutta()` function. A controlled stepper is created by passing the stepper type and the error tolerances to the :code:`odeint::make_controlled` function.

.. code-block:: cpp

    // Absolute and relative tolerance
    double AbsTol = 1e-5, RelTol = 1e-5;

    auto ctrlStepper = odeint::make_controlled<stepper_type>(AbsTol, RelTol);

Adaptivity and more details on the steppers are described in the `boost documentation <https://live.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html>`_.




Performing sensitivity analysis
--------------------------------
We wish to find the sensitivities of the solution at time :math:`t = t_f` w.r.t. the initial conditions :math:`\left[x(t_i),v(t_i)\right]^T` and the damping coefficient :math:`\mu`. The objective functions are the final values of the state vector 

.. math::
    \begin{bmatrix}
        x(t_f)\\
        v(t_f)
    \end{bmatrix}

so :code:`int Nout = 2`.

Forward Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We set the following initial conditions, damping coefficient and times,

.. code-block:: cpp

    // Time interval and time step
    double ti = 0.0 ,tf = 5e-1,dt = 0.001;

    std::vector<double> mu = {1e3};
    std::vector<double> x0 = {2.0, -2.0 / 3.0 + 10.0 / (81.0 * mu[0]) - 292.0 / (2187.0 * mu[0] * mu[0])};

We solve the ODE by performing a forward pass,

.. code-block:: 

    // Create the system function instance
    auto vdp = VanDerPol();

    // Forward pass
    size_t numSteps = runge_kutta(ctrlStepper, vdp, x0, mu, ti, tf, dt, driver);

Notice how the controlled stepper is passed to the :code:`runge_kutta()` function, instead of the stepper itself.

Reverse Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We set the partial derivatives of the solution state vector at time :math:`t=t_f` with respect to the state vector solution of the ODE (:math:`\left[x(t_f),v(t_f)\right]^T`) and with respect to the parameter :math:`\mu`,

.. code-block:: cpp

    auto lambda = std::vector(N, std::vector<double>(N));

    lambda[0][0] = 1.0; // dx(t_f)/d x(t_f) = 1.0
    lambda[0][1] = 0.0; // dx(t_f)/d v(t_f) = 0.0
    lambda[1][0] = 0.0; // dv(t_f)/d x(t_f) = 0.0
    lambda[1][1] = 1.0; // dv(t_f)/d v(t_f) = 1.0

    auto muadj = std::vector(N, std::vector<double>(Npar));
    muadj[0][0] = 0.0; // dx(t_f)/d(mu) = 0
    muadj[1][0] = 0.0; // dv(t_f)/d(mu) = 0

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, muadj);

Additionally, we need to inform the driver of the chosen stepper'same Butcher Tableau,

.. code-block:: 

    stepper_type stepper;

    constructDriverButcherTableau(driver, stepper);

As of now, we still need an instance of the underlying error stepper to inform the Driver of the steppers Butcher Tableau, since I did not find an easy way to access the error stepper from the controlled stepper.

The rest of the code will be the same as in :doc:`harmonicOscillator`.