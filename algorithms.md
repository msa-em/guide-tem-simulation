---
title: Electron Scattering Algorithms
numbering:
  enumerator: 1.%s
math:
  \ii: '{i\mkern1mu}'
  \invFFT: \mathcal{F}^{-1}_{\mathbf{k}\to\mathbf{r}}
  \FFT: \mathcal{F}_{\mathbf{r}\to\mathbf{k}}
  \angstroms: \text{\normalfont\AA}
label : algorithms_page
---


(numerical-solutions-of-the-schrodinger-equation)=
## Numerical Solutions of the Schrödinger Equation

As discussed in [](#physics_bound_systems), the [Schrödinger equation](wiki:Schrödinger_equation) typically cannot be solved analytically in complex systems. Therefore, in order to perform electron scattering simulations, we must calculate numerical solutions of Equation [](#eq:Schrodinger_time) for electron waves. First, we define the {cite:t}`debroglie1925recherches` wavelength of a free electrons (corrected for relativistic effects) as

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

In our simulations, we will assume the {math}`z`-position coordinate of the wavefunction {math}`\psi(\bm{r})` is alone sufficient to describe its propagation in both time and space. Starting from Equaton [](#eq:Schrodinger_time) and assuming steady-state conditions with $V(\bm{r}, t) = V(\bm{r})$, we replace $\psi(\bm{r}, t)$ with $\psi(\bm{r}) e^{-\ii E t / \hbar}$ and separate the time dependence. This gives the time-independent Schrödinger equation:

```{math}
\left[-\frac{\hbar^2}{2m} \nabla^2 + e V(\bm{r}) \right] \psi(\bm{r}) = E \psi(\bm{r}),
```
where $V(\bm{r})$ is the crystal potential, and $E$ is the total energy of the electron. Substituting $\lambda $ and $\sigma$ into the equation and rearranging in terms of the electron’s wavevector $k_0 = 1/\lambda$, we obtain:

```{math}
:label: eq:schrodinger_start
\left[\nabla^2 + 4\pi^2 k_0^2\right] \psi(\bm{r}) = -4\pi^2 \sigma V(\bm{r}) \psi(\bm{r}),
```
This 3D time-independent Schrödinger equation serves as the foundation for deriving numerical electron scattering algorithms. 


(multislice-method)=
## The Multislice Method

By a wide margin, the most common algorithm used for electron scattering simulations is the multislice method, first described by {cite:t}`cowley1957scattering`.
In this method, we make two assumptions:
- The $\partial^2 / \partial z^2$ term in the Laplacian can be neglected, as the wavefunction's variation along the $z$-axis is much slower compared to its variation in the transverse $(x, y)$ directions.
- The longitudinal wavevector $k_0$ is much larger than the contributions from transverse components of the wavefunction, i.e., $k_0 \gg |{\nabla_{xy}}^2|$.

With these assumptions, we can substitute Equations [](#eq:wavelength) and [](#eq:interaction_constant) into Equation [](#eq:schrodinger_start) to obtain {cite:p}`kirkland2020`

```{math}
:label: eq:Shrodinger_electron
\frac{\partial }{\partial z} \psi(\bm{r})
    =
    \frac{\ii \lambda}{4 \pi} {\nabla_{xy}}^2 \psi(\bm{r})
    + 
    \ii \sigma V(\bm{r}) \psi(\bm{r}),
```

where {math}`{\nabla_{xy}}^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2`. 

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
    \psi_0(\bm{r}),
```
where

```{math}
V_{\Delta z}(\bm{r})
    =
    \int_{z_0}^{z_0 + \Delta z} 
    V(\bm{r}) dz,
```

is a thin slice of the potential as described in Equation [](#isolated-atomic-potentials) or [](#dft-potentials).
Unfortunately, even with the above approximations, Equation [](#eq:Shrodinger_simple) cannot be solved in closed form due to the two non-commuting operators. 
Instead, we solve it numerically by using a split-step method, where we alternate between solving each operator independently.
The steps of the multislice method are detailed below.

### Steps of the Multislice Method
#### 1. Atomic Coordinates

We first generate a set of atomic coordinates for our desired sample.
The atomic coordinates are placed in a *simulation cell*, a [rectangular cuboid](#wiki:Rectangular_cuboid) where all edges vectors are $90^\circ$ apart. 
We assume the optic axis of the electron beam is along the $z$ axis. 
For each atom, we define the $\bm{r}=(x,y,z)$ position, the atomic number, a thermal vibration parameter, and sometimes the occupancy.
Ideally the atomic coordinates will be periodic in the $(x,y)$ plane, though this is not always possible.

#### 2. Potential Slices

Next, we calculate the potential $V(\bm{r})$ for the sample. We compute this potential numerically, either using the parameterization approach shown in [](#isolated-atomic-potentials) or using a DFT calculation as described in [](#dft-potentials).
We divide the atomic potentials into *slices*, which are thin sections of the sample in the $(x,y)$ plane. 
Thinner slices will produce more accurate simulations, at the cost of longer computation times.
Typical slice thicknesses for accurate simulations are 1--2 $\rm{\AA}$, equal to roughly the atomic spacing of most solid materials.

When using isolated atomic potentials, we can either take the infinite projected potential which places the full scattering cross-section of each into a single slice, or perform numerical integration of a finite 3D projected potential which allows the potential of each atom to be spread into multiple adjacent slices.

We can also add additional electrostatic or electromagnetic fields to the potential slices. Electrostatic fields can be produced by electric fields across the sample or excess charges or holes, and will produce the same phase shifts as the atomic potentials, described by Equation [](#eq:Shrodinger_electron).
The effect of both extrinsic and intrinsic magnetic fields can be calculated using the [Aharonov–Bohm equation](#wiki:Aharonov–Bohm_effect).

#### 3. Initial Wavefunctions

Next, we define the intitial condition of the electron beam wavefunction $\psi(\bm{r})$, described in [](#CTF_page). In an ideal plane wave TEM or diffraction pattern simulation, we use only a single initial wavefunction. 
We can also include spatial coherence in a [TEM simulation](#tem_sims) by performing a multislice simulation where the initial probe is tilted to a range of incident probe angles, which are then summed incoherently to generate the simulation output.
For a [STEM simulation](#stem_sims), we may need to calculate thousands or even millions of initial conditions for the electron probe, as each unique STEM probe position requires another simulation.

#### 4. Transmission Operator

Following {cite:t}`kirkland2020`, if we assume a slice is infinitesimal thickness, we can set the ${\nabla_{xy}}^2$ term from [](#eq:Shrodinger_simple) to zero and obtain the solution
```{math}
:label: eq:operator_transmission
\psi(\bm{r})
    = 
    \psi_0(\bm{r})
    \exp[\ii \sigma V_{\Delta z}(\bm{r})].
```

We see from this expression that as the electron wavefuncton passes through a given slice, it will pick up a forward phase shift proportional to $V_{\Delta z}(\bm{r})$. This first Born approximation is quite accurate for high accelerating voltages, for small-to-intermediate atomic number species, and for thin slices. However we may require a more accurate expansion and / or numerical slicing of individual atomic potentals when using very low accelerating voltages or for calculating scattering from high atomic number species.

#### 5. Propagation Operator

Next, we need to *propagate* the electron wave from one slice to the next by using the propagation operator in Equation [](#eq:Shrodinger_simple). We assume empty space between slices, setting $V(\bm{r})=0$ in [](#eq:Shrodinger_simple) to get

```{math}
:label: eq:prop01
\psi(\bm{r})
    = 
    \exp \left\{
    \frac{\ii \lambda \Delta z}{4 \pi} {\nabla_{xy}}^2
    \right\}
    \psi_0(\bm{r}).
```
Setting $\Lambda = \lambda \Delta z / 4 \pi$ and Taylor expanding this expression gives
```{math}
:label: eq:prop02
\psi(\bm{r})
    = 
    \left[
      \sum_{m=0}^\infty 
      (\ii \Lambda)^m 
      \frac{\partial^{2m} \psi_0(\bm{r})}{\partial x^{2m}} 
    \right]
    \left[
      \sum_{n=0}^\infty 
      (\ii \Lambda)^n 
      \frac{\partial^{2n} \psi_0(\bm{r})}{\partial y^{2n}} 
    \right].
```
Taking the 2D Fourier transform $\Psi(\bm{k}) = \mathscr{F}_{\bm{r} \rightarrow \bm{k}}\{ \psi(\bm{r}) \}$ of both sides and using the fact that the x and y derivatives are orthogonal, we get
```{math}
:label: eq:prop01
\begin{aligned}
\Psi(\bm{k})

    &= 
    \left[
      \sum_{m=0}^\infty 
      (\ii \Lambda)^m 
      (\ii 2 \pi k_x)^{2m}
    \right]
    \left[
      \sum_{n=0}^\infty 
      (\ii \Lambda)^n 
      (\ii 2 \pi k_y)^{2n}
    \right]
    \Psi_0(\bm{k}) \\
    
    &=
    \left[
      \sum_{m=0}^\infty 
      (-\ii 4 \pi^2 \Lambda {k_x}^2)^m 
    \right]
    \left[
      \sum_{m=0}^\infty 
      (-\ii 4 \pi^2 \Lambda {k_y}^2)^m 
    \right]
    \Psi_0(\bm{k}) \\

    &=
    \left[
      \sum_{m=0}^\infty 
      (-\ii \pi \lambda \Delta z {k_x}^2)^m 
    \right]
    \left[
      \sum_{m=0}^\infty 
      (-\ii \pi \lambda \Delta z {k_y}^2)^m 
    \right]
    \Psi_0(\bm{k}) \\

    &=
    \exp\left(
      -\ii \pi \lambda \Delta z {k_x}^2 
    \right)
    \exp\left(
      -\ii \pi \lambda \Delta z {k_y}^2
    \right)
    \Psi_0(\bm{k}).
\end{aligned}
```
We can now write the final propagation operator by combining ${k_x}^2+{k_y}^2=|\bm{k}|^2$ to get
```{math}
:label: eq:prop
\Psi(\bm{k})
  =
  \exp\left(
    -\ii \pi \lambda \Delta z |\bm{k}|^2 
  \right)
  \Psi_0(\bm{k}).
```

If there are still remaining slices that the electron wave has not passed through, we alternate steps 4 and 5 until the prope wavefunction reaches the output surface of the sample, where it is referred to as the `exit wave`.

#### 6. Transfer Function

After we have calculated the exit wave, we then need to apply the effects of our microscope optics to this wave and reach the detector plane by using a modulation transfer function (MTF). The MTF could be very simple; for example, in either a TEM diffraction simulation or a typical STEM simulation, we assume that the detector is placed at the far field limit and that therefore we only need to Fourier transform the exit wave to reach the detector plane. 

For a TEM imaging simulation, we typically use a contrast transfer function (CTF) for the MTF. The CTF can include aplanatic [optical aberrations](wiki:Optical_aberration) such as defocus, spherical aberration, astigmatism, and higher order coherent wave aberrations. It can also include more complex optical affects such as field distortion, image rotation, or planatic aberrations, where the aberrations vary as a function of position. The CTF equations are described in [](#CTF_page).

#### 7. Detector Functions

Finally, we convert from the complex wavefunction to a real-valued detector measurement. This intensity measurement may be performed in real space for near-field imaging giving $I(\bm{r})$, or in Fourier space for far-field diffraction space measurements giving $I(\bm{k})$. The measured intensity for a pixelated detector is just the magnitude squared of the wavefunction $|\psi(\bm{r})|^2$ or $|\psi(\bm{k})|^2$. To simulated an integrating detector intensity $I_D(\bm{k})$, such as those for BF or DF STEM measurements, we apply a detector function $D(\bm{k})$ to our measured intensity using the expression
```{math}
:label: eq:detector_function
I_D(\bm{R})
    =
    \int_{\bm{k}} 
    |\psi(\bm{R},\bm{k})|^2
    D(\bm{k})
    d\bm{k},
```
where $\bm{R}$ is the position of the STEM probe, and $D(\bm{k})$ is usually an array of zeroes and ones defining the detector shape.

Because we're performing a simulation, we do not need to used a fixed detector geometry. We could instead define variable detectors such as a set of concentric annular ring detectors with a spacing $\Delta k$, using the expression
```{math}
:label: eq:detector_annular_rings
I(\bm{R},n)
    =
    \int_{n \Delta k}^{(n+1) \Delta k} 
    \frac{1}{2 \pi}
    \int_0^{2 \pi} 
    |\psi(\bm{R},\bm{k})|^2
    d\theta dk',
```
where *$\theta$* is the annular coordinate and $k'$ is the radial coordinate for $\bm{k}$-space.

We will use the multislice method in the rest of this article, so you can proceed to [](#wave-aberrations), or read below for information on other simulation methods!

(blochwave-method)=
## Bloch Wave Simulations

While the multislice method is widely used for electron scattering simulations due to its flexibility and scalability, the Bloch wave method provides an alternative approach that is especially efficient for small, periodic structures. 
The Bloch wave formalism leverages the translational symmetry of the crystal lattice to express the electron wavefunction as a sum of periodic Bloch states, significantly reducing computational effort for crystalline samples.
This reduction in computational cost is possible because in a Bloch wave simulation, we only consider a small number of scattering vectors $k$, or equivalently to a small number of scattering angles $\alpha = \lambda k$

The Bloch wave method is particularly advantageous for periodic systems, as it reduces the computational domain to a single unit cell and directly incorporates crystal symmetry. However, it is less suited for non-periodic systems or those with large-scale defects, where the multislice method is more appropriate.
This approach requires careful numerical handling of eigenvalue decomposition and the summation over a sufficiently large number of reciprocal lattice vectors to ensure convergence.

### Steps of the Bloch Wave Method
#### 1. Bloch Wave Expansion

The electron wavefunction $\psi(\bm{r})$ inside a crystal can be expressed as a linear combination of Bloch waves, $b_j(\bm{k}_j, \bm{r})$, which satisfy the periodicity of the crystal potential:

```{math}
\psi(\bm{r}) = \sum_j \alpha_j b_j(\bm{k}_j, \bm{r}),
```
where

```{math}
b_j(\bm{k}_j, \bm{r}) = e^{2\pi i \bm{k}_j \cdot \bm{r}} \sum_{\bm{g}} c_{\bm{g},j} e^{2\pi i \bm{g} \cdot \bm{r}},
```
where
$\bm{k}_j$ are the Bloch wavevectors,
$\bm{g}$ are the reciprocal lattice vectors,
$c_{\bm{g},j}$ are coefficients describing the contribution of each plane wave to the Bloch wave.
This expansion allows us to represent the electron wavefunction as a superposition of states that inherently respect the periodicity of the crystal.

#### 2. Schrödinger Equation

We can rewrite Equation [](eq:schrodinger_start) as:
```{math}
\left[\nabla^2 + 4\pi^2 k_0^2\right] \psi(\bm{r}) = -4\pi^2 \sigma V(\bm{r}) \psi(\bm{r}),
```
where we approximate the second $z$-derivative term for high-energy electrons using $\frac{\partial^2}{\partial z^2} \approx -(2\pi/\lambda)^2$. This approximation assumes the electron wavefunction is dominated by a plane wave propagating along $z$ with wavevector $k_0 = 1 / \lambda$.

To separate the rapidly oscillating component of the wavefunction, we use the substitution $\psi(\bm{r}) = \exp(2\pi i k_0 z) \phi(\bm{r})$, where $\phi(\bm{r})$ represents a slowly varying envelope function. Substituting into the equation yields:

```{math}
\left[-\frac{\hbar^2}{2m} \nabla^2 + eV(\bm{r})\right] \phi(\bm{r}) = E \phi(\bm{r}),
```
where $E = \hbar^2 k_0^2 / 2m$ is the total energy of the electron. This time-independent Schrödinger equation now describes the interaction of the electron wave with the crystal potential $V(\bm{r})$ in all spatial directions.
To account for the periodicity of the crystal lattice, we expand $\phi(\bm{r})$ as a sum of Bloch waves:

```{math}
\phi(\bm{r}) = \sum_{\bm{g}} c_{\bm{g},j} e^{2\pi i \bm{g} \cdot \bm{r}},
```
where $\bm{g}$ are reciprocal lattice vectors, and $c_{\bm{g},j}$ are the coefficients of the expansion. Similarly, the crystal potential is expressed as a Fourier series:

```{math}
V(\bm{r}) = \sum_{\bm{g}} V_{\bm{g}} e^{2\pi i \bm{g} \cdot \bm{r}}.
```
Substituting these expansions into the Schrödinger equation results in a set of coupled equations for the plane wave coefficients $c_{\bm{g},j}$, which form the basis for Bloch wave simulations.

#### 3. Eigenvalue Problem

Inserting the expansions into the Schrödinger equation yields:

```{math}
\sum_{\bm{g}} \left(k_0^2 - |\bm{k}_j + \bm{g}|^2\right) c_{\bm{g},j} e^{2\pi i (\bm{k}_j + \bm{g}) \cdot \bm{r}} 
= -\sum_{\bm{g},\bm{h}} V_{\bm{g} - \bm{h}} C_{\bm{h},j} e^{2\pi i (\bm{k}_j + \bm{g}) \cdot \bm{r}}.
```
By matching coefficients of $\exp^{2\pi i (\bm{k}_j + \bm{g}) \cdot \bm{r}}$, we obtain the eigenvalue equation:

```{math}
\left[2k_0 s_{\bm{g}} - 2\gamma_j k_{0,z}\right] c_{\bm{g},j} + \sum_{\bm{h} \neq \bm{g}} V_{\bm{g} - \bm{h}} C_{\bm{h},j} = 0,
```
where $s_{\bm{g}} = (k_0^2 - |\bm{k}_0 + \bm{g}|^2) / 2k_0$ is the excitation error.
We solve this set of linear equations to find the eigenvalues $2\gamma_j k_{0,z}$ and eigenvectors $c_{\bm{g},j}$, representing the Bloch wave propagation constants and coefficients.

#### 4. Bloch Wave Propagation

The wavefunction $\psi(\bm{r})$ at depth $z$ is expressed as:

```{math}
\psi(\bm{r}) = \sum_{\bm{g}} \psi_{\bm{g}}(z) e^{2\pi i (\bm{k}_0 + \bm{g}) \cdot \bm{r}},
```
where $\psi_{\bm{g}}(z)$ propagates according to:

```{math}
\psi_{\bm{g}}(z) = \sum_j \alpha_j c_{\bm{g},j} e^{2\pi i \gamma_j z}.
```
The propagation constants $\gamma_j$ determine how each Bloch wave evolves through the crystal.

#### 5. Input Wavefunction

At the entrance surface ($z=0$), the wavefunction at the entrance surface of the crystal $\psi_{\bm{g}}(0)$ is matched to the incident wavefunction just outside of the crystal $\psi_0(\bm{r})$,
```{math}
\psi_0(\bm{r}) = \sum_{\bm{g}} \psi_{\bm{g}}(0) e^{2\pi i \bm{g} \cdot \bm{r}},
```
where $\psi_{\bm{g}}(0)$ represents the Fourier components of the input wavefunction. These Fourier components serve as the basis for expanding $\psi_0(\bm{r})$ in terms of Bloch waves.

The expansion of $\psi_0(\bm{r})$ into Bloch waves is achieved by determining the weighting coefficients $\alpha_j$, which describe the contribution of each Bloch wave $b_j(\bm{r})$ to the input wavefunction. By solving:

```{math}
\alpha_j = \sum_{\bm{g}} c_{\bm{g},j}^* \psi_{\bm{g}}(0),
```
we relate the Bloch wave expansion coefficients $\alpha_j$ to the plane wave coefficients $\psi_{\bm{g}}(0)$ of the input wavefunction and the coupling coefficients $c_{\bm{g},j}$, which describe the relationship between the plane wave and Bloch wave bases.

#### 6. Output Wavefunction

At the exit surface ($z = z_{\text{max}}$, corresponding to the crystal thickness $t$), the real-space wavefunction $\psi(\bm{r})$ is expressed using the Bloch wave expansion. The wavefunction is given by:

```{math}
\psi(\bm{r}) = \sum_{j} \alpha_j \left( \sum_{\bm{g}} c_{\bm{g},j} e^{2\pi i \bm{g} \cdot \bm{r}} \right) e^{2\pi i \gamma_j t}.
```
where $\alpha_j$ are the weighting coefficients determined by the input wavefunction, and $c_{\bm{g},j}$ are the Bloch wave coefficients that describe the periodic components of the wavefunction. The plane waves $e^{2\pi i \bm{g} \cdot \bm{r}}$ correspond to reciprocal lattice vectors $\bm{g}$, following the lattice periodicity.
The term $e^{2\pi i \gamma_j t}$ accounts for the propagation of each Bloch wave through the crystal thickness $t$, with $\gamma_j$ representing the eigenvalues that describe the wavevector components along $z$.

Note in particular how the thickness $t$ modulates the contribution of each Bloch wave to the final wavefunction at the exit surface, and how the Bloch wave method can be used to quickly compute output wavefunctions for multiple crystal thicknesses.
