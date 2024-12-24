---
title: Scientific Python Code
numbering:
  enumerator: 1.%s
---

`Python` has inarguably become the leading programming language in the scientific community, continually gaining popularity for the past two decades. 
This is due to `Python` being friendly for beginners, owing to its less strict language structure (dynamic typing, automatic memory management, non-compiled language) and the plethora of tutorials and examples. 
It also benefits from being open-source, with numerous high quality libraries covering a broad range of applications (numeric calculations, image analysis, machine learning, etc.), including a number of packages dedicated to S/TEM. 
The diverse and broad scope of the python ecosystem allows for interoperability between different domains, and file formats as well as easily re-purposing existing implementations algorithms to new domains.

There are a number of `Python` codes relevant to S/TEM simulation and analysis:
- [*ab*TEM](https://github.com/abTEM/abTEM) - all-Python S/TEM image simulation
- [pyPrismatic](https://prism-em.com/tutorial-python/) - image simulation (Python wrapper to C++ package Prismatic)
- [pyMultislice](https://github.com/HamishGBrown/py_multislice) - S/TEM image simulation 
- [py4DSTEM](https://github.com/py4dstem/py4DSTEM) - 4D-STEM analysis
- [liberTEM](https://github.com/LiberTEM/LiberTEM) - 4D-STEM analysis
- [Hyperspy](https://hyperspy.org/) - General framework for S/TEM analysis.
- [pyxem](https://github.com/pyxem/pyxem) - 4D-STEM analysis as part of Hyperspy.
- [Rosetta Scientific Input Output](https://github.com/hyperspy/rosettasciio) - library for microscopy file reading and writing.

The main libraries used in the code examples of this text are:
- NumPy (np) - fast numerical calculations.
- CuPy (cp) - drop-in replacement for NumPy to run on GPUs.
- Numba - just-in-time Python compiler for scientific and array-oriented computing.
- matplotlib - plotting and visualization.
- ipywidgets - making interactive figures.
- Atomic Simulation Environment (ASE) - creating and visualizing atomic structures.

Dynamically interpreted languages including `Python` are particularly attractive for domain experts and scientists trying out new ideas, but the performance of the interpreter can be a barrier to high performance. 
However, by making appropriate use of `Python` open-source numerical libraries including NumPy, CuPy, and Numba, it is possible to write a purely `Python`-based code that performs as well or even better than prior codes based on traditional languages such C++ or Fortran. 
That is indeed what *ab*TEM has been able to achieve, which has made it one of the fastest-growing and popular tools for S/TEM simulations, and our choice for this text.
