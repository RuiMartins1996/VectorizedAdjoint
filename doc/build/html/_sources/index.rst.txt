.. Documentation for VectorizedAdjoint documentation master file, created by
   sphinx-quickstart on Fri Sep 27 16:25:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Tutorial for sensitivity analysis of the Vanderpol equation
============================================================

Introduction
------------
This tutorial demonstrates how to perform sensitivity analysis of the Vanderpol equation using the library provided by this repository.

The Vanderpol equation is a second-order system of nonlinear differential equations 
with two depedent variables :math:`u_1(t)` and :math:`u_2(t)`, given by:

.. math::
   \begin{bmatrix}
	   \dot{u}_1\\
	   \dot{u}_2
   \end{bmatrix}
   =
   \begin{bmatrix}
	   u_2\\
	   \mu\left((1-u_1^2)u_2-u_1\right)
   \end{bmatrix}.

We consider the time domain :math:`[t_i,t_f] = [0,0.5]`, parameter value :math:`\mu=1\times10^3` 
and initial conditions :math:`u_1(0) = 2.0` and :math:`u_2(0) = - 2/3 +10/(81\mu) - 292/(2187\mu^2)`. 


Define the ODE system
--------------------------------
Mathematically, the state of the ODE system is an n-dimensional vector of real or complex numbers (this library only supports real numbers).
In this library, the state of the system is represented by the :code:`std::vector<T>` container type, where :code:`T` is a template parameter.

The ODE system is assumed to be in the form :math:`\dot{x}(t) = f(x,\mu,t)`, where :math:`x` is the state vector, :math:`\mu` is the parameter vector and :math:`t` is the time.
The system is defined by defining a functor that implements the right hand side with the ()-operator with a specific signature:

.. code-block:: cpp

   class VanDerPol
   {
   public:
      VanDerPol() = default;

      template <typename T>
      void operator()(const std::vector<T> &x, std::vector<T> &dxdt, std::vector<T> &mu, const T t) const
      {
         dxdt[0] = x[1];
         dxdt[1] = mu[0] * ((1.0 - x[0] * x[0]) * x[1] - x[0]);
      } 
      
   };

The signature of the ()-operator must always be the same as it is shown in the code above. The first argument is the state vector, the second argument is the rhs vector, the third argument is the parameter vector and the fourth argument is the time.
The ()-operator must be templated to allow an unique implementation to be used by the stepper and by the AADC library. It is also possible to define two versions of the functor, one for :code:`T = double` and another for :code:`T = idouble`, a type that defines an active variable in the AADC library. 
This, however, is uglier code.

Define the Stepper Type
--------------------------------
This library uses the :code:`boost::numeric::odeint` steppers to solve the ODE system. 

.. code-block:: cpp

   typedef odeint::runge_kutta_fehlberg78<std::vector<double>> stepper_type;

In this example we are using an adaptive stepper, the Runge-Kutta-Fehlberg 7(8) stepper, so we need to provide an absolute and relative error tolerance (defined in :code:`const double tol = 1e-5` or by user input in this example).

Performing sensitivity analysis
--------------------------------
If we wish to find the sensitivities of a cost function that depends on the solution of the ODE system and on the parameters, we need to create an instance of the :code:`Driver` class

.. code-block:: cpp
   
   // System size
   int N = 2;
   // Number of parameters
   int Npar = 1;

   // Create a driver object
   Driver driver(N, N, Npar);

This class need to know the size of the state vector, the number of output cost functions and the number of parameters. Here, we are computing the sensitivity matrix of the ODE system, so the output cost functions are the solution of the ODE system (which has a size of 2). 

Performing a forward pass
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   // Forward pass
   size_t numSteps = runge_kutta(
      odeint::make_controlled<stepper_type>(AbsTol, RelTol), vdp, x0, mu, ti, tf, dt, driver);

If we are performing integration with constant step size, the first argument should be an instance of :code:`stepper_type`, where this is an appropriate constant time step stepper.
If we want adaptive step size, we need to pass an instance of a controlled stepper, which can be constructed via :code:`odeint::make_controlled<stepper_type>(AbsTol, RelTol)`, 
where :code:`AbsTol` and :code:`RelTol` are the absolute and relative error tolerances, respectively and :code:`stepper_type` is an error stepper.

Define an Observer
--------------------------------
The Observer construct is a way for the user to log information of the forward and reverse passes. 

.. code-block:: cpp

   class Observer
   {
   public:
      Observer() = default;

      template <typename State, typename Time>
      void operator()(const State &x, Time time)
      {
         std::cout << "t:" << time << "\t";
         std::cout << " x = [";
         for (auto &xi : x) {
            std::cout << xi << " ";
         }
         std::cout << " ]" << std::endl;
      }
   private:
   std::vector<double> times = {};
   };


.. toctree::
   :maxdepth: 2
   :caption: Contents:

