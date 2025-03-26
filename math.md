---
title: Mathematical Concepts
numbering:
  enumerator: 1.%s
---

Here we briefly outline the mathematical concepts required to run STEM simulations. 

(math_numbers)=
## Scalars, Vectors, Tensors, and Arrays

You are probably very familiar with the [real numbers](wiki:Real_number), which are values which measure a continuous quantity. 
These numbers can be integers such as $64$, fractions such as $5/7$, decimals such as $-0.827$ (a subset of fractions), and irrational numbers such as $\sqrt{2}$. 
These values are [scalars](wiki:Scalar_(mathematics)), meaning there is only one number associated with each value, making them one-dimensional (1D).

By contrast, a [vector](wiki:Vector_(mathematics_and_physics)) is a value which cannot be represented by a single number. 
Examples include three-dimensional (3D) position coordinates $(x,y,z)$, magnetic moments, or a set of three Euler angles that describe an arbitrary rotation in 3D space. 
Vectors are always one-dimensional. 
In the same way a vector generalizes the concept of a scalar, a [tensor](wiki:Tensor) generalizes both vectors and scalars by extending a value to any number of dimensions. 
You are also likely implicitly familiar with image [arrays](wiki:Array_(data_structure)), which are 2D tensors for a greyscale image, and 3D for a color image.


(math_complex_numbers)=
## Complex Numbers

To model scattering in electron microscopy, we work with a special number system referred to as [complex numbers](wiki:Complex_number). 
A complex number is a length-2 vector, where the two values are called the **real** and **imaginary** values. 
The real part of each number is the kind of real number as described as above. 
The imaginary part is also a real number, but multiplied by the [imaginary unit](wiki:Imaginary_unit) $i$. 
This unit is a number which satisfies the equation
```{math}
:label: eq:imag
i^2 = -1.
```
Examples of complex numbers include $3+5i$, $4/5-(9/7)i$, $-0.23+4.33i$, and $1+\sqrt{5}i$. [](#fig_complex) shows examples of complex numbers. 
These are enormously powerful, and essential when performing calclations using [quantum mechanics](wiki:Quantum_mechanics). 

```{figure} #app:complex_numbers
:name: fig_complex
:placeholder: ./figures/complex_numbers.png
**Complex numbers.** Every complex value consists of a real and an imaginary component, which can be used to compute the magnitude and phase of each value.
```

In addition to the real and imaginary parts, we can define two other useful quantities. 
The magnitude of a complex number $z=a+ib$ is defined by its Euclidean length
```{math}
:label: eq:mag
|z| = \sqrt{a^2 + b^2}.
```
This quantity is also sometimes referred to as the magnitude or the absolute value. 
The second quantity we define is the complex argument, defined as the angle $\phi$ from the real axis to the complex vector. 
Note that we usually want the signed angle, meaning we use which quadrant the complex value $(a,b)$ is located in to determine the signed angle using the expression
```{math}
:label: eq:arg
\phi = \arctan2(b, a).
```
In programming, functions which measure this value are often named $\arg(z)$ or $\rm{angle}(z)$. 
One of the most powerful identities using complex numbers is [Euler's formula](wiki:Euler%27s_formula)
```{math}
:label: eq:euler
\exp(i z) = \cos(z) + i \sin(z).
```
This expression links complex numbers to the trigonometric functions. 
This relationship foreshadows a key upcoming concept, namely that the trigonometric oscillations describe waves, where the physical concept of the complex angle will be mapped onto relative phase shifts of electron waves. 
We will therefore usually refer to the complex argument as the **phase** of a given complex value.

As we will see, even though computers internally store complex numbers using real and imaginary components as $z=(a,b)$, it is often more intuitive for us to represent these numbers using their magnitude and phase, i.e. as $|z| \exp{(i \phi)}$.

We also sometimes take the [conjugate](wiki:Complex_conjugate) of complex numbers, which means that we reverse the sign of the imaginary part. 
This operation is defined for a complex number $z=a+i b$ as
```{math}
:label: eq:complex_conjugate
z^* = a - i b
```
or for $z=|z| \exp{(i \phi)}$ as 
```{math}
:label: eq:complex_conjugate_2
z^* = |z|\exp{(-i \phi)}.
```

(math_fourier_transform)=
## 1D Fourier Transform

The [Fourier transform](wiki:Fourier_transform) is an integral transformation that decomposes a function into its constituent spatial frequencies. 
This means that we write a function as a series of orthogonal [sine waves](wiki:Sine_wave) with different frequencies / periodicities. 
A set of functions are [orthogonal](wiki:Orthogonal_functions) if we can multiply any pair of them and integrate over their domain to get a result of zero. 
The key property conveyed by this orthogonality is that the set of sine waves output by the Fourier Transform form a complete basis, and thus can represent any arbitrary function.  

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

(math_discrete_fourier_transform)=
## 1D Discrete Fourier Transform

The Fourier transform is very useful for many mathematical analyses. 
However, we cannot use it in this form for this article, because unlike continuous variables such as position $x$ and frequency $\omega$, our experimental images and simulated electron waves are sampled on discrete grids. 
We instead use the [discrete Fourier transform](wiki:Discrete_Fourier_transform), which transforms a finite-length function that has been sampled on a uniformly-spaced grid $f(r)$ into its Fourier transform $F(k)$ using the expression
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
where $N$ is the number of samples for both functions, and $f_n$ is the $n$'th sample of $f(x)$. 
The inverse DTF is defined as
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
where $F_k$ is the $k$'th sample of $F(k)$. 
Both $f(x)$ and $F(k)$ can be complex-valued. 
If either sequence is purely real-valued, the other will be conjugate symmetric, which means that if reverse the sequence and take the complex conjugate, we will recover the same sequence.

[](#fig:fft_1d) shows an interactive version of the 1D discrete Fourier transform. 
Try switching between the different preset functions. 
You can also change the amplitude of $F(k)$ by clicking and moving the points up and down, or the phase by moving left and right. 
It can be difficult at first to gain intution for how changing the values of a given spatial frequency $F_k$ will change the output, but working in Fourier space will become second nature to you when you routinely perform STEM simulations.

```{figure} #app:fft_1d
:name: fig:fft_1d
:placeholder: ./static/fft_1d.png
**Interactive discrete 1D Fourier transform.**
```


(math_fast_fourier_transform)=
## 2D Fourier Transforms

When simulating transmission electron microscopy images, diffraction patterns, and other wavefunctions, we are typically working with two-dimensional (2D) arrays.
Thus we need to expand our previous definition to work in two dimensions. 
However, because this transformation is orthogonal over different Cartesian coordinates, we can compute a 2D discrete Fourier transform by simply sequentially multiplying the two transformations, giving
```{math}
:label: eq:dft_2d_xy
\begin{aligned}
  F(k_x, k_y) 
  &= \mathscr{F}_{x \rightarrow k_x}\{
    \mathscr{F}_{y \rightarrow k_y}\{
      f(x,y)
    \} 
  \} 
  \\
  &= \sum_{m=0}^{M-1} \sum_{n=0}^{N-1} 
    f_{nm}
    \exp\left(
      -2 i \pi k_x m / M
    \right) 
    \exp\left(
      -2 i \pi k_y n / N
    \right), 
\end{aligned}
```
where $f_{nm}$ is the pixel at indices $(m,n)$ for the function $f(x,y)$, which has $(M,N)$ total number of pixels in the $x$ and $y$ directions respectively. 
We can simplify this notation by combining the two sequential operations into one operator, giving the forward and inverse transforms
```{math}
:label: eq:dft_2d
\begin{aligned}
  F(\bm{k}) 
  &= \mathscr{F}_{r \rightarrow k}\{
      f(\bm{r})
  \} 
  \\
  f(\bm{r}) 
  &= \mathscr{F}_{k \rightarrow r}\{
      F(\bm{k})
  \},
\end{aligned}
```
where $\bm{r} = (x,y)$ and $\bm{k} = (k_x,k_y)$ are the 2D coordinate systems for real and diffraction space respectively. 
As in the 1D case, in practice we use the [fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) algorithm to compute these transforms with small computational overhead. 
[](#fig:fft_2d) shows an interactive version of the forward and and inverse 2D FFT. 


```{figure} #app:fft_2d
:name: fig:fft_2d
:placeholder: ./static/fft_2d.png
**Interactive discrete 2D Fourier transform.** Left panel shows diffraction space, right shows real space upsampled by 4x. 
Both images have been FFT-shifted, with origins shown as a red `+` symbol. 
Amplitude is shown as the pixel values, while phase is shown as a periodic colorwheel.
```

On the left panel of [](fig:fft_2d), you can click to swap individual pixels between 0 and 1 in Fourier space, while observing the response on the real space function. 
On the right panel, you can click and drag to shift the function by a distance $\Delta \bm{r}$, which is applied by using the Fourier shift theorem as
```{math}
f(\bm{r} - \Delta \bm{r}) =
  \mathscr{F}_{k \rightarrow r}\{
    \mathscr{F}_{r \rightarrow k}\{
      f(\bm{r})
    \}
    \exp( -2i \pi \bm{k} \, \cdot \Delta \bm{r} )  
  \},  
``` 
where $\exp( -2i \pi \bm{k} \, \cdot \Delta \bm{r} )$ is the expression corresponding to a plane wave. 
If you try dragging the real space function in [](#fig:fft_2d) a short distance, you will see the characteristic coloring of a plane wave, a periodically repeating rainbow function.

A deep understanding of the 2D Fourier transform is essential to both simulate and understand TEM experiments. 
Therefore in [](#tab:fft_pairs) we provide a table of Fourier transform pairs which are important for TEM. 
Note that for simplicity we have omitted some details such as some of the numerical prefactors in front of each function.


(math_fourier_pairs)=
## Fourier transform pairs

:::{list-table} 2D Fourier transform pairs important in TEM simulation and analysis. Note that these transform pairs assume that both $\bm{r}$ and $\bm{k}$ are 2D, i.e. $\bm{r}=(x,y)$ and $\bm{k}=(k_x,k_y)$.
:label: tab:fft_pairs
:header-rows: 2

* - 
  - Diffraction space
  - Real space
  - 

* - Name
  - Function
  - Function
  - Name

* - Delta Dirac 
  - $\delta(\bm{k})$
  - 1
  - plane wave along optic axis

* - Shifted Dirac function
  - $\delta(\bm{k} - \bm{k}_0)$
  - $\exp(-2i \pi \, \bm{r} \cdot \bm{k}_0)$
  - plane wave tilted to $\bm{k}_0$

* - diffraction spot pair
  - $\delta(\bm{k} - \bm{k}_0) + \delta(\bm{k} + \bm{k}_0)$
  - $\cos(-2 \pi \, \bm{r} \cdot \bm{k}_0)$
  - cosine wave


* - plane wave mult.
  - $F(\bm{k}) \exp(-2i \pi \bm{k} \cdot \Delta \bm{r})$
  - $f(\bm{r} - \Delta \bm{r})$
  - Shifted function

* - circular aperture
  - $u(|\bm{k}| - k_{\rm{max}})$
  - $\frac{J_1(2 \pi k_{\rm{max}} |\bm{r}| )}{\pi k_{\rm{max}} |r|}$
  - Airy disk function

* - defocus operator
  - $\exp(-i\pi \lambda |\bm{k}|^2 C_1)$
  - $\exp
  \left(
    i \frac{\pi |\bm{r}|^2}{\lambda C_1}
    \right)$
  - Fresnel quadratic phase factor

* - Gaussian envelope
  - $\exp
  \left(
    -\frac{|\bm{k}|^2}{2 \sigma^2}
  \right)$
  - $\exp
  \left(
    -2 \pi^2 \sigma^2 |\bm{r}|^2
  \right)$
  - Gaussian envelope

* - exponential decay
  - $\exp(-\lambda |\bm{k}|)$
  - $\frac{\lambda^3}{(\lambda^2 + |\bm{r}|^2)^{3/2}}$
  - 3D Lorentzian function

* - Lorentzian function
  - $\frac{\gamma^2}{\gamma^2 + |\bm{k}|^2}$
  - $K_0(2 \pi \gamma |\bm{r}|)$
  - Modified Bessel func., $2^{\rm{nd}}$ kind
:::

Additional transform pairs can be found in the [Wikipedia table of important Fourier transforms](wiki:Fourier_transform#Tables_of_important_Fourier_transforms).

<!-- <table>
  <tr>
    <td>One</td>
    <td>Two</td>
  </tr>
  <tr>
    <td colspan="2">Three</td>
  </tr>
</table> -->
