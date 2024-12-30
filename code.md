---
title: Scientific Python Code
numbering:
  enumerator: 1.%s
---

## Overview

`Python` has unarguably become the leading programming language in the scientific community, continually gaining popularity for the past two decades. 
This is due to `Python` being friendly for beginners, owing to its less strict language structure (dynamic typing, automatic memory management, non-compiled language) and the plethora of tutorials and examples. 
It also benefits from being open-source, with numerous high-quality libraries covering a broad range of applications (numeric calculations, image analysis, machine learning, etc.), including a number of packages dedicated to STEM. 
The diverse and broad scope of the python ecosystem allows for interoperability between different domains and file formats as well as easy re-purposing of existing algorithms to new domains.

There are a number of `Python` codes relevant to STEM simulation and analysis:
- [*ab*TEM](https://github.com/abTEM/abTEM) - all-Python STEM image simulation
- [pyPrismatic](https://prism-em.com/tutorial-python/) - STEM image simulation (Python wrapper to C++ package Prismatic)
- [pyMultislice](https://github.com/HamishGBrown/py_multislice) - STEM image simulation 
- [py4DSTEM](https://github.com/py4dstem/py4DSTEM) - 4D-STEM analysis
- [liberTEM](https://github.com/LiberTEM/LiberTEM) - 4D-STEM analysis
- [Hyperspy](https://hyperspy.org/) - general framework for STEM analysis
- [pyxem](https://github.com/pyxem/pyxem) - 4D-STEM analysis as part of Hyperspy
- [Rosetta Scientific Input Output](https://github.com/hyperspy/rosettasciio) - library for microscopy file reading and writing.

The main libraries used in the code examples of this text are:
- NumPy (np) - fast numerical calculations.
- CuPy (cp) - drop-in replacement for NumPy to run on GPUs.
- Numba (nb) - just-in-time Python compiler for scientific and array-oriented computing.
- matplotlib - plotting and visualization.
- ipywidgets - making interactive figures.
- Atomic Simulation Environment (ASE) - creating and visualizing atomic structures.

Dynamically interpreted languages including `Python` are particularly attractive for domain experts and scientists trying out new ideas, but the performance of the interpreter can be a barrier to high performance. 
However, by making appropriate use of `Python` open-source numerical libraries including NumPy, CuPy, and Numba, it is possible to write a purely `Python`-based code that performs as well or even better than prior codes based on traditional languages such C++ or Fortran. 
That is indeed what *ab*TEM has been able to achieve, which has made it one of the fastest-growing and popular tools for STEM simulations, and our choice for this text.


## Simple Python Examples

Most python scripts begin with importing necessary libraries. 
A *library* is a collection of pre-written code that provides tools, functions, and methods to perform specific tasks, such as mathematical operations, data visualization, or file manipulation, without requiring the user to write everything from scratch.
We can import libraries and alias them (meaning we substitute as a short form name to reduce typing) using this syntax:

```python
import numpy as np
import matplotlib.pyplot as plt
```

`Numpy` is a very useful library for computing mathematical functions. 
For example if we define a coordinate system $r=(x,y)$ where both $x$ and $y$ have a range -5 to +5 in steps of 0.1, and we want to evaluate the function $f(\vec{r}) = \exp(-|r|^2/(2 \sigma^2)$, we could type

```python
sigma = 2.0
x = np.arange(-5,5,0.1)
y = np.arange(-5,5,0.1)

ya,xa = np.meshgrid(y,x)

f = np.exp(
  -(xa**2 + ya**2) / (2*sigma**2)
)
```

where we have set $\sigma=2$. 
We have used the function `np.arange` to calculate a range from [-5,+5) in steps of 0.1, which does not include the final point at +5, instead ending at $5 - 0.1 = 4.9$. 
We have then used the function `np.meshgrid` to expand our 1D $x$ and $y$ vectors into 2D arrays, where the $x$ coordinates are repeated along the $y$-axis and vice versa. 
Finally we compute the exponential function using `np.exp`. 
Note the indentation inside the exponential function, which we use to improve readability.

One of the key features of 'Numpy' arrays is that they can be *broadcasted*, which means smaller arrays are automatically expanded in size during operations to match the shape of larger arrays, without creating unnecessary memory copies or requiring explicit loops.
We can therefore reshape $x$ and $y$ to be a row and column vector, allowing us to calculate the same expression as above without needing to call `np.meshgrid`:

```python
sigma = 2.0
x = np.arange(-5,5,0.1)[:,None]
y = np.arange(-5,5,0.1)[None,:]

f = np.exp(
  -(x**2 + y**2) / (2*sigma**2)
)
```

Here we modify the $x$ vector, by using square brackets $[,]$ to access indices of the $x$ array. 
Specifically we use `[:,None]`, where the first `:` means "all elements" and places these along the 1st dimension (array rows), and the `None` (or alternatively `np.newaxis`) expands the array to become 2D, by adding a new 2nd dimension (array columns). 
This changes the $x$ array shape `x.shape` from $(200,)$ to $(200,1)$. 
The similar indexing applied to $y$ changes its array shape `y.shape` from $(200,)$ to $(1,200)$. 
Finally, when using the addition operator, `Numpy` *broadcasts* these two arrays from shapes $(200,1)$ to $(1,200)$ to $(200,200)$, which gives the compatible array shapes and allows them to be added together.

We might also want to plot the results of our calculation. 
We can use the `Matplotlib` plotting library to plot this image, for example with the code which produces the result shown in [](#fig:gaussian_plotting_example):

```python
plt.figure(figsize=(4,3))
plt.imshow(f, 
  extent=(-5, 5, -5, 5), 
  origin='lower', 
  cmap='viridis',
)
plt.colorbar(label='Amplitude')
plt.title("2D Gaussian Function")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
```

```{figure} figures/gaussian_plotting_example.jpg
:name: fig:gaussian_plotting_example
Example `python` plotting of a 2D Gaussian distribution.
```

Note that here we have moved the default origin for arrays from the upper left corner down to the lower left corner, using the argument `origin='lower'`. 
This change of origin is a common source of errors for `Python` newcomers, as in mathematics generally (and `Numpy` specifically) we typically place the origin in the upper left corner, while most plotting functions generally (and `Matplotlib` specifically) place it in the lower left. 
We should remember however that the directions are only display and plotting conventions; internally the entries in our 2D `Numpy` arrays are always self-consistent when being accessed using the `[index_0,index_1]` convention.

## Python Notebooks

'Python' notebooks, such as those popularized by 'Jupyter', are a versatile tool for scientific computation and exploration. 
They allow users to combine code, text, mathematics (using LaTeX), and visualizations in a single document. 
This format is particularly advantageous for iterative workflows where the ability to test new ideas, visualize intermediate results, and document findings is critical. 
Each section of code, called a *cell*, can be executed independently, which makes debugging and refining specific parts of a larger workflow straightforward. 
The downside of this organizational scheme is that we must carefully check that our code works correctly when all cells are executed sequentially.

For tools like *ab*TEM, notebooks provide an excellent platform to run simulations, tweak their input parameters, and visualize the results.
For example, users can modify the defocus or aperture size in a simulated STEM image and then immediately see how the results change.
This approach not only makes learning the software more intuitive but also facilitates efficient experimentation, enabling users to focus on understanding the science.
The ability to share notebooks further enhances collaboration and reproducibility, as researchers can publish notebooks alongside their papers, allowing readers to verify the results and modify them for their own needs.

In this article, we use computable Jupyter notebooks alongside the text to ensure reproducibility.
You can inspect all the source code under the `Notebooks` directory on the table of contents.
In addition, we use the `ipywidgets` and `ipympl` libraries to make the code *interactive* using sliders and cursor callback functions.


