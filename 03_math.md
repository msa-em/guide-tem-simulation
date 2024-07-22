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
i^2 = -1.
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
|z| = \sqrt{a^2 + b^2}.
```
This quantity is also sometimes referred to as the magnitude or the absolute value. The second quantity we define is the complex argument, defined as the angle $\phi$ from the real axis to the complex vector. Note that we usually want the signed angle, meaning we use which quadrant the complex value $(a,b)$ is located in to determine the signed angle using the expression
```{math}
:label: eq:arg
\phi = \arctan2(b, a).
```
In programming, functions which measure this value are often named $\arg(z)$ or $\rm{angle}(z)$. One of the most powerful identities using complex numbers is [Euler's formula](wiki:Euler%27s_formula)
```{math}
:label: eq:euler
\exp(i z) = \cos(z) + i \sin(z).
```
This expression links complex numbers to the trigonometric functions. This relationship  foreshadows a key upcoming concept where the trigonometric oscillations describe waves, where the physical concept of the complex angle will be mapped onto relative phase shifts of electron waves. We will therefore usually refer to the complex argument as the **phase** of a given complex value.

As we will see, even though computers internally store complex numbers using real and imaginary components as $z=(a,b)$, it is often more intiutive for us to represent these numbers using their magnitude and phase, i.e. as $|z| \exp{(i \phi)}$.

We also sometimes take the [complex conjugate](wiki:Complex_conjugate) of complex numbers, which means that we reverse the sign of the imaginary part. This operation is defined for a complex number $z=a+i b$ as
```{math}
:label: eq:complex_conjugate
z^* = a - i b
```
or for $z=|z| \exp{(i \phi)}$ as 
```{math}
:label: eq:complex_conjugate_2
z^* = \exp{(-i \phi)}.
```



(fourier_transform)=
## 1D Fourier Transform

The [Fourier Transform](wiki:Fourier_transform) is an integral transform that decomposes a function into its constituent spatial frequencies. This means that we write a function as a series of orthogonal [sine waves](wiki:Sine_wave) with different frequencies / periodicities. A series of functions are [orthogonal](wiki:Orthogonal_functions) if we can multiply any pair of them and integrate over their domain to get a result of zero. The key property conveyed by this orthogonality is that the set of sine waves output by the Fourier Transform form a complete basis, and thus can represent any arbitrary function.  

We can map a function $f(x)$ over the $x$ coordinates to its Fourier transform $F(\omega)$ over the $\omega$ coordinates using the expression
```{math}
:label: eq:ft_forward
F(\omega) = \int_{-\infty}^{\infty}
  f(x) \exp\left(
    -i \omega x
  \right) dx,
```
and we can recover the original function with the inverse transform
```{math}
:label: eq:ft_inverse
f(x) = \frac{1}{2 \pi} \int_{-\infty}^{\infty}
  F(\omega) \exp\left(
    -i \omega x
  \right) d\omega.
```


(discrete_fourier_transform)=
## 1D Discrete Fourier Transform

The Fourier transform is very useful for many mathematical analyses. However, we cannot use it in this form for this paper, because unlike the continuous variables $x$ and $\omega$, our experimental images and simulated electron waves are sampled on discrete grids. We instead use the [Discrete Fourier transform](wiki:Discrete_Fourier_transform) (DFT), which transforms a finite-length function which has been samplied on a uniformly-spaced grid $f(r)$ into its Fourier transform $F(k)$ using the expression
```{math}
:label: eq:dft_forward
\begin{aligned}
  F(k) 
  &= \mathscr{F}_{r \rightarrow k}\{ f(r) \} 
  \\
  &= \sum_{n=0}^{N-1} 
    f_n \exp\left(
      -2 i \pi k n / N
    \right), 
\end{aligned}
```
where $N$ is the number of samples for both functions, and $f_n$ is the $n$'th sample of $f(x)$. The inverse DTF is defined as
```{math}
:label: eq:dft_inverse
\begin{aligned}
  f(r) 
  &= \mathscr{F}_{k \rightarrow r}\{ F(k) \} 
  \\
  &= \frac{1}{N} \sum_{k=0}^{N-1} 
    F_k \exp\left(
      2 i \pi k n / N
    \right),
\end{aligned}
```
where $F_k$ is the $k$'th sample of $F(k)$. Both $f(x)$ and $F(k)$ can be complex-valued. If either sequence is purely real-valued, the other will be conjugate symmetric, which means that if reverse the sequence and take the complex conjugate, we will recover the same sequence.

[](#fig:fft_1d) shows an interactive version of the 1D DFT. Try switching between the different preset functions. You can also change the amplitude of $F(k)$ by clicking and moving the points up and down, or the phase by moving left and right. It can be difficult at first to gain intution for how changing the values of a given spatial frequency $F_k$ will change the output, but working in Fourier space will become second nature to you if you routinely perform S/TEM simulations.


```{figure} #app:fft_1d
:name: fig:fft_1d
:placeholder: ./static/fft_1d.png
**Discrete 1D Fourier Transform.**
```


(fast_fourier_transform)=
## 2D Fast Fourier Transform


