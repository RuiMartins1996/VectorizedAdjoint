Tutorial for sensitivity analysis of an harmonic oscillator
============================================================

Introduction
------------
This tutorial demonstrates how to perform sensitivity analysis of 
a damped harmonic oscillator using the library provided by this repository.

The model can be described by the following system of ODEs:

.. math::
   \begin{bmatrix}
	   \dot{x}\\
	   \dot{v}
   \end{bmatrix}
   =
   \begin{bmatrix}
	   v\\
	   -kx - \mu v
   \end{bmatrix},

where :math:`x` is the displacement relative to the equilibrium position, :math:`v` is the velocity, :math:`k` is the spring constant and :math:`\mu` is the damping coefficient.

Define the ODE system
--------------------------------
Mathematically, the state of the ODE system is an :math:`n`-dimensional vector of real or complex numbers (this library only supports real numbers).
In this library, the state of the system is represented by the :code:`std::vector<T>` container type, where :code:`T` is a template parameter.

The ODE system is assumed to be in the form :math:`\dot{\mathbf{r}}(t) = f(\mathbf{r},\alpha,t)`, where :math:`\mathbf{r}` is the state vector, :math:`\alpha` is the parameter vector and :math:`t` is the time. We will assume that :math:`k = 1`, so: 

.. math::
    \mathbf{r} =
    \begin{bmatrix}
        x\\
        v
    \end{bmatrix}
    \qquad
    \alpha =
    \begin{bmatrix}
        \mu
    \end{bmatrix}

The ODE system is defined by defining a functor that implements the right hand side :math:`f` with the ()-operator with a specific signature:

.. code-block:: cpp

    class HarmonicOscillator
    {
    private:
        double k = 1.0;

    public:
        HarmonicOscillator() = default;

        template <typename T>
        void operator()(const std::vector<T> &r,std::vector<T> &drdt,const std::vector<T> &mu,const T t) const
        {
            drdt[0] = r[1];
            drdt[1] = - k * r[0] - mu[0] * r[1];
        }
    };

Notice how parameters that we don't wish to differentiate with respect to can be defined as class attributes (:code:`double k = 1.0`), while the parameters we wish to differentiate with respect to are passed to the ()-operator as :code:`std::vector<T> &mu`.

The signature of the ()-operator must always be the same as it is shown in the code above. The ()-operator must be templated to allow the definition of the system function only once, since a version of the system function with :code:`T = double` is needed by the stepper and a version with :code:`T = idouble` is need by the AADC library to record the rhs with automatic differentiation. It is also possible to define two versions of the functor, one for :code:`T = double` and another for :code:`T = idouble` and pass each version to the correct place. This, however, is uglier code.

Define the Stepper Type
--------------------------------
This library uses the :code:`boost::numeric::odeint` steppers to solve the ODE system. 

.. code-block:: cpp

    typedef odeint::runge_kutta4<std::vector<double>> stepper_type;

In this example we are using a constant step stepper, the classic 4th-ordr Runge-Kutta method.

Performing sensitivity analysis
--------------------------------
Assume that we wish to find the sensitivities of the mechanical energy of the harmonic oscillator at :math:`t = t_f` with respect to the initial conditions :math:`\left[x(t_i),v(t_i)\right]^T` and damping coefficient :math:`\mu`. The mechanical energy, which is the objective function, is given by:

.. math::
    E(t_f) = \frac{1}{2}kx(t_f)^2 + \frac{1}{2}mv(t_f)^2

where :math:`m = 1.0`. 

We need to create an instance of the :code:`Driver` class,

.. code-block:: cpp
   
   // System size
   int Nin = 2;
   // Number of cost functions 
   int Nout = 1;
   // Number of parameters
   int Npar = 1;

   Driver driver(Nin, Nout, Npar);

which is the construct that stores the trajectory of the ODE system. 

Forward Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We need to set the initial conditions, initial and final time, time step and the value of the parameter :math:`\mu`:

.. code-block:: cpp

    // Time interval and time step
    double ti = 0.0, tf = 10.0, dt = 0.01;

    std::vector<double> mu = {0.15};
    std::vector<double> r0 = {0.0, 1.0};

For this choice of damping coefficient, the harmonic oscillator is underdamped. We then solve the ODE by doing a forward pass,

.. code-block:: 

    // Create a stepper instance
    stepper_type stepper;
    // Create the system function instance
    auto hm = HarmonicOscillator();

    // Forward pass
    size_t numSteps = runge_kutta(stepper, hm, r0, mu, ti, tf, dt, driver);

and the solution is stored in :code:`r0`.

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



