.. Documentation for VectorizedAdjoint documentation master file, created by
   sphinx-quickstart on Fri Sep 27 16:25:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Main Page 
===============

This is the main page. Below are links to documentation regarding the examples organized by name. They should be followed in order.

The tutorials will demonstrate how to perform sensitivity analysis of several models using the library provided by this repository.

We will be considering problems where we have an ODE system given in the form 

.. math:: 
   \dot{\mathbf{x}}(t) = f(\mathbf{x},\alpha,t)`
   
where :math:`\mathbf{r}` is the state vector, :math:`\alpha` is the parameter vector and :math:`t` is the time.

All the library function and classes are contained within the :code:`vectorizedadjoint` namespace.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   harmonicOscillator
   vanderpol
   generalizedLotkaVolterra
