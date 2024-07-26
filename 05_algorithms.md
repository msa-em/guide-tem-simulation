---
title: Electron Scattering Algorithms
numbering:
  enumerator: 1.%s
math:
  \ii: '{i\mkern1mu}'
  \invFFT: \mathcal{F}^{-1}_{\mathbf{k}\to\mathbf{r}}
  \FFT: \mathcal{F}_{\mathbf{r}\to\mathbf{k}}
  \angstroms: \text{\normalfont\AA}
---


(numerical-solutions-of-the-schrodinger-equation)=
### Numerical Solutions of the Schrödinger Equation

As discussed in [](#physics_page), the [Schrödinger equation](wiki:Schrödinger_equation) typically cannot be solved analytically in complex systems. Therefore, in order to perform electron scattering simulations, we must calculate numerical solutions of [](#eq:Schrodinger_time) for electron waves. First, we define the {cite:t}`debroglie1925recherches` wavelength of a free electrons (corrected for relativistic effects) as

```{math}
:label: eq:wavelength

\lambda = \frac{h \, c}{\sqrt{e \, E_0 (2 \, m \, c^2 + e \, E_0)}},
```

where {math}`h` is the Plank constant, {math}`c` is the speed of light, {math}`e` is the electron charge, and {math}`E_0` is the accelerating voltage applied to the electron. Using SI units for [](#eq:wavelength) will give the wavelength in units of meters. In practice we typically use length units of  for all calculations, and therefore multiple this result by {math}`10^{10}`.

Next, we define the electron-potential interaction constant as (the numerical values of these constants can be found in  [](#app:constants))

```{math}
:label: eq:interaction_constant

\sigma = \frac{2 \pi \, m \, e \, \lambda}{h^2}.
```

In our simulations, we will assume the {math}`z`-position coordinate of the wavefunction {math}`\psi(\bm{r})` is alone sufficient to describe its propagation in both time and space, and therefore drop the {math}`t` coordinate. Substituting [](#eq:wavelength) and [](#eq:interaction_constant) into [](#eq:Schrodinger_time), we obtain {cite:p}`kirkland2020`

```{math}
:label: eq:Shrodinger_electron
\frac{\partial }{\partial z} \psi(\bm{r})
    =
    \frac{\ii \lambda}{4 \pi} {\nabla_{xy}}^2 \psi(\bm{r})
    + 
    \ii \sigma V(\bm{r}) \psi(\bm{r}),
```

where {math}`{\nabla_{xy}}^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2`. To perform electron scattering simulations, we numerically solve Equation [](#eq:Shrodinger_electron) using one of the methods described below.



(multislice-method)=
### The Multislice Method

By a wide margin, the most common algorithm used for electron scattering simulations is the multislice method, first described by {cite:t}`cowley1957scattering`.
Equation [](#eq:Shrodinger_electron) shows the overall numerical recipe we will use; when the wavefunction {math}`\psi_0(\bm{r})` is at position {math}`z_0`, we will evaluate the operators on the right hand side over a distance {math}`\Delta z` to calculate the new wavefunction {math}`\psi(\bm{r})` at position {math}`z_0 + \Delta z`. {cite:t}`kirkland2020` gives the formal operator solution to [](#eq:Shrodinger_electron) as

```{math}
:label: eq:Shrodinger_solution

\psi(\bm{r})
    = 
    \exp \left\{
    \int_{z_0}^{z_0 + \Delta z} 
    \left[
        \frac{\ii \lambda}{4 \pi} {\nabla_{xy}}^2
        + 
        \ii \sigma V(\bm{r})
    \right] dz
    \right\}
    \psi_0(\bm{r})
```

Assuming {math}`\Delta z` is small, [](#eq:Shrodinger_solution) can be simplified to

```{math}
:label: eq:Shrodinger_simple
\psi(\bm{r})
    = 
    \exp\left[
        \frac{\ii \lambda}{4 \pi} \Delta z {\nabla_{xy}}^2
        + 
        \ii \sigma V_{\Delta z}(\bm{r})
    \right]
    \psi_0(\bm{r})
```

where

```{math}
V_{\Delta z}(\bm{r})
    =
    \int_{z_0}^{z_0 + \Delta z} 
    V(\bm{r}) dz,
```

is a thin slice of the potential as described in [](#isolated-atomic-potentials) or [](#dft-potentials).
Unfortunately, even with the above approximations, [](#eq:Shrodinger_simple) cannot be solved in closed form due to the two non-commuting operators. 
Instead, we solve it numerically by using a split-step method, where we alternate between solving each operator independently.
The steps of the multislice method are detailed below.


#### 1 - Atomic Coordinates

We first generate a set of atomic coordinates for our desired sample.
The atomic coordinates are placed in a *simulation cell*, a [rectangular cuboid](#wiki:Rectangular_cuboid) where all edges vectors are $90^\circ$ apart. 
We assume the optic axis of the electron beam is along the $z$ axis. 
For each atom, we define the $\bm{r}=(x,y,z)$ position, the atomic number, a thermal vibration parameter, and sometimes the occupancy.
Ideally the atomic coordinates will be periodic in the $(x,y)$ plane, though this is not always possible.


#### 2 - Potential Slices

Next, we calculate the potential $V(\bm{r})$ for the sample. We compute this potential numerically, either using the parameterization approach shown in [](#isolated-atomic-potentials), or using a DFT calculation as described in [](#dft-potentials).
We divide the atomic potentials into *slices*, which are thin sections of the sample in the $(x,y)$ plane. 
Thinner slices will produce more accurate simulations, at the cost of longer computation times.
Typical slice thicknesses for accurate simulations are 1-2 $\rm{\AA}$, equal to roughly the atomic spacing of most solid materials.

When using isolated atomic potentials, we can either take the infinite projected potential which places the full scattering cross-section of each into a single slice, or perform numerical integration of a finite 3D projected potential which allows the potential of each atom to be spread into multiple adjacent slices.

We can also add additional electrostatic or electromagnetic fields to the potential slices. Electrostatic fields can be produced by electric fields across the sample or excess charges or holes, and will produce the same phase shifts as the atomic potentials, described by Equation [](#eq:Shrodinger_electron).
The effect of both extrinsic and intrinsic magnetic fields can be calculated using the [Aharonov–Bohm equation](#wiki:Aharonov–Bohm_effect).


#### 3 - Initial Wavefunctions

Next, we define the intitial condition of the electron beam wavefunction $\psi(\bm{r})$, described in Sections `insert sections`. In an ideal plane wave TEM or diffraction pattern simulation, we use only a single initial wavefunction. 
We can also include spatial coherence in a [TEM simulation](#tem_sims) by performing a multislice simulation where the initial probe is tilted to a range of incident probe angles, which are then summed incoherently to generate the simulation output.
For a [STEM simulation](#stem_sims), we may need to calculate thousands or even millions of initial conditions for the electron probe, as each unique STEM probe position requires another simulation.


#### 4 - Transmisson Operator

Following {cite:t}`kirkland2020`, if we assume a slice is infinitesimal thickness, we can set the ${\nabla_{xy}}^2$ term from [](#eq:Shrodinger_simple) to zero and obtain the solution
```{math}
:label: eq:operator_transmission
\psi(\bm{r})
    = 
    \psi_0(\bm{r})
    \exp[\ii \sigma V_{\Delta z}(\bm{r})].
```

We see from this expression that as the electron wavefuncton passes through a given slice, it will pick up a forward phase shift proportional to $V_{\Delta z}(\bm{r})$. This first Born approximation is quite accurate for high accelerating voltages, for small-to-intermediate atomic number species, and for thin slices. However we may require a more accurate expansion and / or numerical slicing of individual atomic potentals when using very low accelerating voltages or for calculating scattering from high atomic number species.


#### 5 - Propagation Operator

Next, we need to *propagate* the electron wave from one slice to the next by using the propagation operator in Equation [](#eq:Shrodinger_simple).




If there are still remaining slices that the electron wave has not passed through, we alternate steps 3 and 4 until the prope wavefunction reaches the output surface of the sample, where it is referred to as the `exit wave`.

#### 6 - Transfer Function

After we have calculated the exit wave, we then need to apply the effects of our microscope optics to this wave and reach the detector plane by using a microscope transfer function (MTF). The MTF could be very simple; for example in either a TEM diffraction simulation or a typical STEM simulation we assume the detector is placed at the far field limit, and therefore only need to Fourier transform the exit wave to reach the detector plane. 

For a TEM imaging simulation, we typically use a contrast transfer function (CTF) for the MTF. The CTF can include aplanatic [optical aberrations](wiki:Optical_aberration) such as defocus, spherical aberration, astigmatism, and higher order coherent wave aberrations. It can also include more complex optical affects such as field distortion, image rotation, or planatic aberrations, where the aberrations vary as a function of position. The CTF equations are described in `add section link`.

#### 7 - Detector Functions

Finally .





Unfortunately, even with the above approximations, [](#eq:Shrodinger_simple)




 is to solve the *transmission* operator by calculating the interaction of the electron beam inside a thin slice

is to find a subset of the atoms 

 take a thin *slice* of the sample along the beam direction by selecting a subset of atoms 






This algorithm is known as the multislice method because 




(bloch-wave-method)=
### The Bloch Wave Method

text


(prism-method)=
### PRISM

text


(other-method)=
### Other Methods

text



<!-- ```{figure} #app:fft_1d

```

```{figure} #app:probe_size

```
 -->