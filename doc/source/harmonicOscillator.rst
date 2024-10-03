Damped Harmonic Oscillator
============================================================

Objective
------------
To find how the mechanical energy at :math:`t = t_f` varies with respect to the initial conditions and the damping coefficient.

Model
------------
A damped harmonic oscillator can be described by the following second-order linear ordinary differential equation:

.. math:: 
    m\ddot{x} + \mu\dot{x} + kx = 0

where :math:`x` is the displacement relative to the equilibrium position, :math:`m` is the mass, :math:`\mu` is the damping coefficient and :math:`k` is the spring constant.

In this example, we will take :math:`k = 1 \text{ }(N/m)` as a constant parameter and :math:`m = 1 \text{ }(kg)` as the mass. 

This equation can be cast into the following system of ODEs:

.. math::
   \begin{bmatrix}
	   \dot{x}\\
	   \dot{v}
   \end{bmatrix}
   =
    f(\begin{bmatrix}x\\v\end{bmatrix},\mu,t)
   =
   \begin{bmatrix}
	   v\\
	   -kx - \mu v
   \end{bmatrix},

where :math:`v` is the velocity. 

The mechanical energy of the oscillator at a given time :math:`t` is given by:

.. math::
    E(t) = \frac{1}{2}kx(t)^2 + \frac{1}{2}mv(t)^2

Data Types 
--------------------------------
Mathematically, the state of a (real valued) ODE system is a :math:`\mathbb{R}^N` vector. In this library, toth state of the system and the parameter vector are represented by the :code:`std::vector<double>` container type. In this example, the state vector and the parameter vector are:

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

time is represented by a scalar of type :code:`double`.

Definition of the ODE system 
--------------------------------
A functor that implements the right hand side of the ODE system, :math:`f(\cdot,\cdot,\cdot)`, with the ()-operator has to be defined:

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

The library expects the ()-operator of the user defined functor to have a specific signature, which is shown in the code above, indepedently of the model being considered.

Notice how parameters that we don't wish to differentiate w.r.t. can be defined as class attributes (:code:`double k = 1.0`), while the parameters we wish to differentiate w.r.t. are passed to the ()-operator in :code:`std::vector<T> &mu`.

The ()-operator should be templated to allow the user to define a single version of the system function, since a version with :code:`T = double` is needed by the stepper and a version with :code:`T = idouble` is need by the AADC library to record the rhs with automatic differentiation. It is also possible to define two versions of the functor, one for :code:`T = double` and another for :code:`T = idouble` and pass each version to the correct place. This, however, is uglier code.

Define the Stepper Type
--------------------------------
This library uses the :code:`boost::numeric::odeint` steppers to solve the ODE system. 

.. code-block:: cpp

    using namespace boost::numeric::odeint;
    typedef odeint::runge_kutta4<std::vector<double>> stepper_type;

In this example we are using a constant time step stepper, the classic 4th-order Runge-Kutta method. More details about steppers can be found in the `boost documentation <https://live.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html>`_.





Performing sensitivity analysis
--------------------------------
We wish to find the sensitivities of the mechanical energy at :math:`t = t_f` w.r.t. the initial conditions :math:`x(t_i),v(t_i)` and damping coefficient :math:`\mu`. 

We need to create an instance of the :code:`Driver` class,

.. code-block:: cpp
   
   // System size
   int Nin = 2;
   // Number of cost functions 
   int Nout = 1;
   // Number of parameters
   int Npar = 1;

   Driver driver(Nin, Nout, Npar);

which is the construct that stores the trajectory of the ODE system and some system characteristics.

Forward Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We need to set the initial conditions, initial and final time, time step and the value of the parameter :math:`\mu`:

.. code-block:: cpp

    // Time interval and time step
    double ti = 0.0, tf = 10.0, dt = 0.01;

    std::vector<double> mu = {0.15};
    std::vector<double> r0 = {0.0, 1.0};

We solve the ODE by performing a forward pass,

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
We set the partial derivatives of the cost function :math:`È(t_f)` with respect to the state vector solution of the ODE (:math:`\left[x(t_f),v(t_f)\right]^T`) and with respect to the parameter :math:`\mu`,

.. code-block:: cpp

    auto lambda = std::vector(Nout, std::vector<double>(Nin));
    lambda[0][0] = 1.0 * r0[0]; // dE(t_f)dx(t_f) = k * x(t_f)
    lambda[0][1] = 1.0 * r0[1]; // dE(t_f)dv(t_f) = m * v(t_f)

    auto muadj = std::vector(Nout, std::vector<double>(Npar));
    muadj[0][0] = 0.0; // dE(t_f)d\mu = 0.0

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, muadj);

Additionally, we need to inform the driver of the chosen stepper's Butcher Tableau,

.. code-block:: 

    constructDriverButcherTableau(driver, stepper);

and to record the rhs function with automatic differentiation (currently, an hand-written alternative has not been implemented),

.. code-block:: 

    recordDriverRHSFunction(driver, hm);


We perform a reverse pass by doing:

.. code-block::

    // Reverse pass to obtain the adjoints of the cost functions
    adjointSolve(driver, mu);

The sensitivities are stored in :code:`lambda` and :code:`muadj`.



