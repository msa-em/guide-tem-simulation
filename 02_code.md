---
title: Scientific Python Code
numbering:
  enumerator: 1.%s
---

Python has inarguably become the leading programming language in the scientific community, continually gaining popularity for the past two decades. This is due to Python being friendly for beginners, owing to its less strict language structure (dynamic typing, automatic memory management, non-compiled language) and the plethora of tutorials and examples. It also benefits from being free (this text would be less accessible if it required a licence), with numerous high quality libraries covering a broad range of applications (numeric calculations, image analysis, machine learning, ...), including a number of packages dedicated to S/TEM. The diverse and broad scope of the python ecosystem allows for interoperability between different domains, and file formats as well as easily re-purposing existing implementations algorithms to new domains.

There are a number of libraries dedicated to S/TEM:
- *ab*TEM - All-Python S/TEM image simulation
- pyPrismatic - Image simulation (Python wrapper to C++ package Prismatic)
- pyMultislice - Image simulation 
- py4DSTEM - 4D-STEM analysis
- libreTEM - 4D-STEM analysis
- pyxem - 4D-STEM
- Hyperspy - General S/TEM analysis 
- ...

The main libraries used here are:
- NumPy (np) - fast numerical calculations
- CuPy (cp) - drop-in replacement for NumPy to run on GPUs
- Numba - just-in-time Python compiler for scientific and array-oriented computing
- matplotlib - all things plotting 
- ipywdigets - making interactive figures
- ASE (Atomic Simulation Environment) - creating and visualizing atomic structures   
- ...

Dynamically interpreted languages including Python are particularly attractive for domain experts and scientists trying out new ideas, but the performance of the interpreter can be a barrier to high performance. However, by making appropriate use of Python open-source numerical libraries including NumPy, CuPy, and Numba, it is possible to write a purely Python-based code that performs as well or even better than prior codes based on traditional languages such C++ or Fortran. That is indeed what *ab*TEM has been able to achieve, which has made it one of the fastest-growing and popular tools for S/TEM simulations, and our choice for this text.