---
title: Mathematic Concepts
numbering:
  enumerator: 1.%s
---

Here we briefly outline the mathematical concepts required to run S/STEM simulations. 

(numbers)=
## Scalars, Vectors, Tensors, and Arrays

You are probably very familiar with the [real numbers](wiki:Real_number), which are values which measure a continuous quantity. These numbers can be integers such as $64$, fractions such as $5/7$, decimals such as $-0.827$ (a subset of fractions), and irrational numbers such as $\sqrt{2}$. These values are [scalars](wiki:Scalar_(mathematics)), meaning there is only one number associated with each value.

By contrast, a [vector](wiki:Vector_(mathematics_and_physics)) is a value which cannot be represented by a single number. Examples include 3D position coordinates $(x,y,z)$, magnetic moments, or a set of 3 Euler angles which describe an arbitrary rotation in 3D space. Vectors are always one-dimensional. In the same way a vector generalizes the concept of a scalar, a [tensor](wiki:Tensor) generalizes both vectors and scalars by extending a value to any number of dimensions. You are also likely familiar with image [arrays](wiki:Array_(data_structure)), which are 2D tensors for a greyscale image, and 3D for a color image.



(complex_numbers)=
## Complex Numbers

To model scattering in electron microscopy, we work with a special number system referred to as [complex numbers](wiki:Complex_number). A complex number is a length-2 vector, where the two values are called the **real** and **imaginary** values. The real part of each number is exactly the same kind of real number as described as above. The imaginary part is also a real number, but multiplied by the [imaginary unit](wiki:Imaginary_unit) $i$. This unit is a number which satisfies the equation
```{math}
:label: eq:imag
\begin{align} 
  i^2 = -1
\end{align}.
```
Examples of complex numbers include $3+5i$, $4/5-(9/7)i$, $-0.23+4.33i$, and $1+\sqrt{5}i$. [](#fig_complex) shows examples of complex numbers. Complex numbers are enormously powerful, and essential when performing calclations using [quantum mechanics](wiki:Quantum_mechanics). 


```{figure} 
:name: fig_complex
:placeholder: ./static/complex_numbers.png
**Complex numbers.** Every complex value consists of a real and an imaginary component, which can be used to compute the magnitude and phase of each value.
```

In addition to the real and imaginary parts, we can define two other useful quantities. The magnitude of a complex number $z=a+ib$ is defined by its Euclidean length
```{math}
:label: eq:mag
\begin{align} 
  |z| = \sqrt{a^2 + b^2}.
\end{align}
```
This quantity is also sometimes referred to as the magnitude or the absolute value. The second quantity we define is the complex argument, defined as the angle $\phi$ from the real axis to the complex vector. Note that we usually want the signed angle, meaning we use which quadrant the complex value $(a,b)$ is located in to determine the signed angle using the expression
```{math}
:label: eq:arg
\begin{align} 
  \phi = \arctan2(b, a).
\end{align}
```
In programming, functions which measure this value are often named $\arg(z)$ or $\rm{angle}(z)$. 

One of the most powerful identities using complex numbers is [Euler's formula](wiki:Euler%27s_formula)
```{math}
:label: eq:euler

\begin{align} 
  \exp(i z) = \cos(z) + i \sin(z)
\end{align}.
```
This expression links complex numbers to the trigonometric functions. This relationship  foreshadows a key upcoming concept where the trigonometric oscillations describe waves, where the physical concept of the complex angle will be mapped onto relative phase shifts of electron waves. We will therefore usually refer to the complex argument as the **phase** of a given complex value.

As we will see, even though computers internally store complex numbers using real and imaginary components as $z=(a,b)$, it is often more intiutive for us to represent these numbers using their magnitude and phase, i.e. as $|z| \exp{(i \phi)}$.




(fourier_transform)=
## The Fourier Transform

text
