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

As discussed in [](#physics_bound_systems), the [Schrödinger equation](wiki:Schrödinger_equation) typically cannot be solved analytically in complex systems. Therefore, in order to perform electron scattering simulations, we must calculate numerical solutions of Equation [](#eq:Schrodinger_time) for the fast electrons in the microscope. In this section we derive the equation solved by the multislice method, and state each approximation as we make it. More detailed derivations can be found in {cite:t}`kirkland2020`.

First, we define the {cite:t}`debroglie1925recherches` wavelength of the incident electrons, corrected for relativistic effects {cite:p}`humphreys_scattering_1979`, as

```{math}
:label: eq:wavelength

\lambda = \frac{h \, c}{\sqrt{e \, E_0 (2 \, m \, c^2 + e \, E_0)}},
```

where {math}`h` is the Planck constant, {math}`c` is the speed of light, {math}`e` is the elementary charge, {math}`m` is the electron rest mass, and {math}`E_0` is the accelerating voltage. Using SI units for [](#eq:wavelength) will give the wavelength in units of meters. In practice we typically use length units of Å for all calculations, and therefore multiply this result by {math}`10^{10}`. For example, the wavelength of 200 kV electrons is $\lambda \approx 0.0251$ Å.

Next, we define the electron-potential interaction constant as (the numerical values of these constants can be found in [](#app:constants))

```{math}
:label: eq:interaction_constant

\sigma = \frac{2 \pi \, \gamma \, m \, e \, \lambda}{h^2},
```

where $\gamma = 1 + e E_0 / (m c^2)$ is the relativistic Lorentz factor, which converts the electron rest mass $m$ into its relativistic mass $\gamma m$. The interaction constant has units of radians per volt per Å, and converts the projected potential of a specimen into a phase shift of the electron wave.

We look for steady-state solutions, where the electron energy is fixed by the accelerating voltage and the potential does not change in time, i.e. $V(\bm{r}, t) = V(\bm{r})$. Starting from Equation [](#eq:Schrodinger_time), we replace $\psi(\bm{r}, t)$ with $\psi(\bm{r}) e^{-\ii E t / \hbar}$ and separate the time dependence. This gives the time-independent Schrödinger equation:

```{math}
:label: eq:schrodinger_tise
\left[-\frac{\hbar^2}{2 \gamma m} \nabla^2 - e V(\bm{r}) \right] \psi(\bm{r}) = E \psi(\bm{r}),
```
where $V(\bm{r})$ is the electrostatic potential of the specimen in volts, $-e V(\bm{r})$ is the potential energy of the electron (negative because the electron is attracted to the positive potential near the atomic nuclei), and $E$ is the kinetic energy of the incident electron far from the specimen. Substituting $\lambda$ and $\sigma$ into the equation and rearranging in terms of the magnitude of the electron's wavevector $k_0 = 1/\lambda$, we obtain:

```{math}
:label: eq:schrodinger_start
\left[\nabla^2 + 4\pi^2 k_0^2\right] \psi(\bm{r}) = -\frac{4\pi \sigma}{\lambda} V(\bm{r}) \psi(\bm{r}).
```
This 3D time-independent Schrödinger equation serves as the foundation for deriving numerical electron scattering algorithms. 


(multislice-method)=
## The Multislice Method

By a wide margin, the most common algorithm used for electron scattering simulations is the multislice method, first described by {cite:t}`cowley1957scattering`.
For fast electrons, the wavefunction is dominated by the incident plane wave travelling along the optic axis $z$. We therefore write the wavefunction as the product

```{math}
:label: eq:envelope
\psi(\bm{r}) = \phi(\bm{r}) \exp\left(2 \pi \ii z / \lambda\right),
```

where $\phi(\bm{r})$ is an envelope function that varies slowly along $z$ compared to the wavelength $\lambda$ {cite:p}`kirkland2020`. Substituting [](#eq:envelope) into [](#eq:schrodinger_start), the $4 \pi^2 k_0^2$ terms cancel, leaving

```{math}
{\nabla_{xy}}^2 \phi(\bm{r}) 
+ \frac{\partial^2 \phi(\bm{r})}{\partial z^2} 
+ \frac{4 \pi \ii}{\lambda} \frac{\partial \phi(\bm{r})}{\partial z}
= -\frac{4\pi \sigma}{\lambda} V(\bm{r}) \phi(\bm{r}),
```

where {math}`{\nabla_{xy}}^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2` is the Laplacian in the transverse $(x,y)$ plane. We now make the paraxial approximation: because $\phi(\bm{r})$ varies slowly along $z$, the term $\partial^2 \phi / \partial z^2$ is much smaller than $(4\pi/\lambda) \, \partial \phi / \partial z$, and can be neglected. This approximation also neglects backscattering, which is very weak for the high accelerating voltages used in S/TEM. Rearranging the remaining terms gives {cite:p}`kirkland2020`

```{math}
:label: eq:Shrodinger_electron
\frac{\partial }{\partial z} \psi(\bm{r})
    =
    \frac{\ii \lambda}{4 \pi} {\nabla_{xy}}^2 \psi(\bm{r})
    + 
    \ii \sigma V(\bm{r}) \psi(\bm{r}),
```

where we have renamed the envelope function $\phi(\bm{r})$ back to $\psi(\bm{r})$. The factor $\exp(2 \pi \ii z / \lambda)$ is identical at every transverse position $(x,y)$, and so it does not change any measured intensity. From this point on we therefore drop it, and $\psi(\bm{r})$ refers to the slowly-varying envelope of the electron wave.

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
    \psi_0(\bm{r}).
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

is the projected potential of a thin slice of the specimen, computed from the isolated atomic potentials described in [](#physics_isolated_atomic_potentials) or the DFT potentials described in [](#physics_dft_potentials).
Unfortunately, even with the above approximations, Equation [](#eq:Shrodinger_simple) cannot be evaluated directly, because the two operators in the exponent do not commute. 
Instead, we solve it numerically by using the split-step method, where we alternate between applying each operator separately, i.e. we approximate $\exp(A + B) \approx \exp(A)\exp(B)$ for the propagation operator $A$ and the transmission operator $B$, with an error that decreases as the slice thickness $\Delta z$ is reduced.
The same split-step approach is widely used to solve other wave equations, including nonlinear Schrödinger equations {cite:p}`taha_analytical_1984`.
The steps of the multislice method are detailed below.

### Steps of the Multislice Method
#### 1. Atomic Coordinates

We first generate a set of atomic coordinates for our desired sample.
The atomic coordinates are placed in a *simulation cell*, a [rectangular cuboid](wiki:Rectangular_cuboid) where all edge vectors are $90^\circ$ apart. 
We assume the optic axis of the electron beam is along the $z$ axis. 
For each atom, we define the $\bm{r}=(x,y,z)$ position, the atomic number, a thermal vibration parameter (see [](#frozen-phonons)), and sometimes the occupancy.
Ideally the atomic coordinates will be periodic in the $(x,y)$ plane, though this is not always possible.

#### 2. Potential Slices

Next, we calculate the potential $V(\bm{r})$ for the sample. We compute this potential numerically, either using the parameterization approach shown in [](#physics_isolated_atomic_potentials) or using a DFT calculation as described in [](#physics_dft_potentials).
We divide the atomic potentials into *slices*, which are thin sections of the sample in the $(x,y)$ plane. 
Thinner slices will produce more accurate simulations, at the cost of longer computation times.
Typical slice thicknesses for accurate simulations are 1–2 Å, equal to roughly the atomic spacing of most solid materials.

When using isolated atomic potentials, we can either take the infinite projected potential which places the full scattering cross-section of each atom into a single slice, or perform numerical integration of a finite 3D projected potential which allows the potential of each atom to be spread into multiple adjacent slices.

We can also add additional electrostatic or magnetic fields to the potential slices. Electrostatic fields can be produced by electric fields across the sample or excess charges or holes, and will produce the same phase shifts as the atomic potentials, described by Equation [](#eq:Shrodinger_electron).
The phase shifts produced by both extrinsic and intrinsic magnetic fields can be calculated from the magnetic vector potential, using the [Aharonov–Bohm effect](wiki:Aharonov–Bohm_effect) {cite:p}`aharonov_significance_1959`.

#### 3. Initial Wavefunctions

Next, we define the initial condition of the electron beam wavefunction $\psi(\bm{r})$, described in [](#CTF_page). In an ideal plane wave TEM or diffraction pattern simulation, we use only a single initial wavefunction. 
We can also include spatial coherence in a [TEM simulation](#tem_sims) by performing a multislice simulation where the initial probe is tilted to a range of incident probe angles, which are then summed incoherently to generate the simulation output.
An incoherent sum adds the intensities $|\psi(\bm{r})|^2$ of the individual simulations, while a coherent sum adds the complex wavefunctions $\psi(\bm{r})$ before the squared magnitude is taken. 
Waves summed coherently can interfere with each other, while waves summed incoherently cannot.
For a [STEM simulation](#stem_sims), we may need to calculate thousands or even millions of initial conditions for the electron probe, as each unique STEM probe position requires another simulation.

#### 4. Transmission Operator

Following {cite:t}`kirkland2020`, if we assume a slice has infinitesimal thickness, we can set the ${\nabla_{xy}}^2$ term from [](#eq:Shrodinger_simple) to zero and obtain the solution
```{math}
:label: eq:operator_transmission
\psi(\bm{r})
    = 
    \psi_0(\bm{r})
    \exp[\ii \sigma V_{\Delta z}(\bm{r})].
```

We see from this expression that as the electron wavefunction passes through a given slice, it will pick up a forward phase shift proportional to the projected potential $V_{\Delta z}(\bm{r})$. 
[](#eq:operator_transmission) is known as the projection, or phase grating, approximation, because it treats all of the scattering within a slice as occurring in a single infinitesimally thin plane. 
Note that this is not the same as the first Born approximation described in [](#physics_born_approximation): the transmission operator includes scattering to all orders of the potential within the slice, but it neglects the propagation of the wave inside the slice. 
Multiple scattering through the full thickness of the specimen is then included by alternating the transmission and propagation steps.
The projection approximation is quite accurate for high accelerating voltages, for small-to-intermediate atomic number species, and for thin slices. However we may require thinner slices and / or numerical slicing of individual atomic potentials when using very low accelerating voltages or for calculating scattering from high atomic number species.

#### 5. Propagation Operator

Next, we need to *propagate* the electron wave from one slice to the next. We assume empty space between slices, setting $V(\bm{r})=0$ in [](#eq:Shrodinger_electron) to get the free-space paraxial wave equation

```{math}
:label: eq:prop01
\frac{\partial }{\partial z} \psi(\bm{r})
    =
    \frac{\ii \lambda}{4 \pi} {\nabla_{xy}}^2 \psi(\bm{r}).
```

This equation is easiest to solve in Fourier space. We take the 2D Fourier transform $\Psi(\bm{k}, z) = \mathscr{F}_{\bm{r} \rightarrow \bm{k}}\{ \psi(\bm{r}) \}$ of both sides, with the sign convention of [](#math_fast_fourier_transform). Each derivative $\partial / \partial x$ becomes a multiplication by $2 \pi \ii k_x$ (and likewise for $y$), and so the transverse Laplacian becomes a multiplication by $(2 \pi \ii k_x)^2 + (2 \pi \ii k_y)^2 = -4 \pi^2 |\bm{k}|^2$. This gives

```{math}
:label: eq:prop02
\frac{\partial }{\partial z} \Psi(\bm{k}, z)
    =
    -\ii \pi \lambda |\bm{k}|^2 \, \Psi(\bm{k}, z),
```

which is an ordinary differential equation in $z$ for each spatial frequency $\bm{k}$. Its solution over a distance $\Delta z$ is

```{math}
:label: eq:prop
\Psi(\bm{k}, z_0 + \Delta z)
  =
  \exp\left(
    -\ii \pi \lambda \Delta z |\bm{k}|^2 
  \right)
  \Psi(\bm{k}, z_0).
```

[](#eq:prop) is the Fresnel free-space propagator written in Fourier space: each spatial frequency picks up a phase shift proportional to $|\bm{k}|^2$, which is the same form as the defocus operator in [](#tab:fft_pairs). In real space, applying this propagator is equivalent to convolving the wavefunction with the Fresnel kernel, which is proportional to $\exp\left[\ii \pi |\bm{r}|^2 / (\lambda \Delta z)\right]$. In practice we always apply the propagator in Fourier space, which requires only two FFTs and an element-wise multiplication.

If there are still remaining slices that the electron wave has not passed through, we alternate steps 4 and 5 until the wavefunction reaches the exit surface of the sample, where it is referred to as the exit wave.

#### 6. Transfer Function

After we have calculated the exit wave, we then need to apply the effects of the microscope optics after the specimen, and propagate the wave to the detector plane. This step can be very simple; for example, in either a TEM diffraction simulation or a typical STEM simulation, we assume that the detector is placed at the far field limit and that therefore we only need to Fourier transform the exit wave to reach the detector plane. 

For a TEM imaging simulation, we apply the contrast transfer function (CTF) of the objective lens to the exit wave. The CTF can include aplanatic [optical aberrations](wiki:Optical_aberration) such as defocus, spherical aberration, astigmatism, and higher order coherent wave aberrations. It can also include more complex optical effects such as field distortion, image rotation, or off-axial aberrations, where the aberrations vary as a function of position in the field of view. The CTF equations are described in [](#CTF_page).

Real detectors also reduce the contrast at high spatial frequencies, due to their finite pixel size and the spreading of the signal from each electron into neighboring pixels. This loss of contrast is described by the modulation transfer function (MTF) of the detector, which we can apply by multiplying the Fourier transform of the measured intensity, after it has been computed in the next step.

#### 7. Detector Functions

Finally, we convert from the complex wavefunction to a real-valued detector measurement. This intensity measurement may be performed in real space for near-field imaging giving $I(\bm{r})$, or in Fourier space for far-field diffraction space measurements giving $I(\bm{k})$. The measured intensity for a pixelated detector is just the magnitude squared of the wavefunction $|\psi(\bm{r})|^2$ or $|\psi(\bm{k})|^2$. To simulate an integrating detector intensity $I_D(\bm{k})$, such as those for BF or DF STEM measurements, we apply a detector function $D(\bm{k})$ to our measured intensity using the expression
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

Because we're performing a simulation, we do not need to use a fixed detector geometry. We could instead define variable detectors such as a set of concentric annular ring detectors with a spacing $\Delta k$, using the expression
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

(inelastic-scattering)=
### Inelastic Scattering

The multislice method as described above simulates elastic scattering, where the fast electron does not lose energy to the specimen. 
Real specimens also scatter electrons inelastically, by exciting phonons (atomic vibrations), plasmons (collective oscillations of the valence electrons), and single-electron transitions such as inner-shell ionization. 
Each of these processes can be included in simulations, with different levels of approximation {cite:p}`allen_modelling_2015`:

- **Phonons.** Thermal diffuse scattering from phonons is usually included using the frozen phonon method, described in [](#frozen-phonons), which averages the intensities of many elastic simulations with different random atomic displacements {cite:p}`loane_thermal_1991`. The energy lost to phonons is very small (tens to hundreds of meV), and so frozen phonon simulations treat phonon scattering as quasi-elastic {cite:p}`van_dyck_frozen_2009,forbes_thermal_2011`.
- **Absorptive potentials.** When we only need the elastic signal, for example in an energy-filtered image, the electrons removed by inelastic scattering can be approximated by adding a small imaginary component to the potential, which reduces the amplitude of the elastic wave {cite:p}`weickenmeier_computation_1991`.
- **Plasmons.** Plasmon scattering can be added to the multislice method using a Monte Carlo approach, where each simulated electron can randomly excite a plasmon in each slice, changing its energy and deflecting it by a small angle {cite:p}`mendis_inelastic_2019`.
- **Inner-shell ionization.** Core-loss signals such as EELS and EDX maps can be simulated by computing the transition potential for each ionization event, which describes the scattering of the fast electron while it excites a core electron to an unoccupied state. The inelastically scattered wave is then propagated through the rest of the specimen using the multislice method {cite:p}`dwyer_multislice_2005`. This requires one additional multislice calculation for each transition and each scattering depth, which can be accelerated substantially using scattering matrix methods {cite:p}`brown_linear-scaling_2019`.

In this article we only consider elastic and thermal diffuse scattering, but the same simulation codes can often be extended to these inelastic signals; for example, *ab*TEM can compute core-loss transitions from DFT calculations.

We will use the multislice method in the rest of this article, so you can proceed to [](#wave-aberrations), or read below for information on other simulation methods!

(blochwave-method)=
## Bloch Wave Simulations

While the multislice method is widely used for electron scattering simulations due to its flexibility and scalability, the Bloch wave method provides an alternative approach that is especially efficient for small, periodic structures. 
The Bloch wave method was introduced by {cite:t}`bethe_theorie_1928`, and detailed descriptions can be found in {cite:t}`humphreys_scattering_1979` and {cite:t}`zuo_dynamical_2017`.
The Bloch wave formalism leverages the translational symmetry of the crystal lattice to express the electron wavefunction as a sum of periodic Bloch states, significantly reducing computational effort for crystalline samples.
This reduction in computational cost is possible because in a Bloch wave simulation, we only consider a small number of scattering vectors $\bm{g}$, or equivalently a small number of scattering angles $\alpha = \lambda |\bm{g}|$.

The Bloch wave method is particularly advantageous for periodic systems, as it reduces the computational domain to a single unit cell and directly incorporates crystal symmetry. However, it is less suited for non-periodic systems or those with large-scale defects, where the multislice method is more appropriate.
This approach requires careful numerical handling of eigenvalue decomposition and the summation over a sufficiently large number of reciprocal lattice vectors to ensure convergence.

### Steps of the Bloch Wave Method
#### 1. Bloch Wave Expansion

The electron wavefunction $\psi(\bm{r})$ inside a crystal can be expressed as a linear combination of Bloch waves, $b_j(\bm{k}_j, \bm{r})$, which satisfy the periodicity of the crystal potential:

```{math}
:label: eq:bloch_sum
\psi(\bm{r}) = \sum_j \alpha_j b_j(\bm{k}_j, \bm{r}),
```
where

```{math}
:label: eq:bloch_wave
b_j(\bm{k}_j, \bm{r}) = e^{2\pi i \bm{k}_j \cdot \bm{r}} \sum_{\bm{g}} c_{\bm{g},j} e^{2\pi i \bm{g} \cdot \bm{r}},
```
where
$\bm{k}_j$ are the Bloch wavevectors,
$\bm{g}$ are the reciprocal lattice vectors,
$c_{\bm{g},j}$ are coefficients describing the contribution of each plane wave to the Bloch wave, and $\alpha_j$ are the excitation coefficients that describe how strongly each Bloch wave is excited by the incident beam.
This expansion allows us to represent the electron wavefunction as a superposition of states that inherently respect the periodicity of the crystal.
Below we show how to compute $\bm{k}_j$, $c_{\bm{g},j}$, and $\alpha_j$.

#### 2. Schrödinger Equation for a Crystal

Inside a crystal, the potential is periodic, and so it can be written as a Fourier series over the reciprocal lattice vectors $\bm{g}$,

```{math}
V(\bm{r}) = \sum_{\bm{g}} V_{\bm{g}} e^{2\pi i \bm{g} \cdot \bm{r}}.
```
We also define the scaled Fourier coefficients $U_{\bm{g}} = \sigma V_{\bm{g}} / (\pi \lambda)$, so that the right hand side of [](#eq:schrodinger_start) becomes $-4 \pi^2 \sum_{\bm{g}} U_{\bm{g}} e^{2\pi i \bm{g} \cdot \bm{r}} \psi(\bm{r})$.
Following the same idea as the multislice derivation, we write the wavefunction as a sum of plane waves, one for each diffracted beam $\bm{g}$, whose amplitudes $\psi_{\bm{g}}(z)$ vary slowly with the depth $z$:

```{math}
:label: eq:bloch_beams
\psi(\bm{r}) = \sum_{\bm{g}} \psi_{\bm{g}}(z) e^{2\pi i (\bm{k}_0 + \bm{g}) \cdot \bm{r}},
```
where $\bm{k}_0$ is the incident wavevector with magnitude $k_0 = 1/\lambda$. For simplicity, we only include the reciprocal lattice vectors $\bm{g}$ which are perpendicular to the beam direction (the zero order Laue zone).

#### 3. Eigenvalue Problem

We substitute [](#eq:bloch_beams) into [](#eq:schrodinger_start), multiply both sides by $e^{-2\pi i (\bm{k}_0 + \bm{h}) \cdot \bm{r}}$, and integrate over the unit cell. 
Because plane waves with different reciprocal lattice vectors are orthogonal over the unit cell, only the terms with $\bm{g} = \bm{h}$ survive this integration, which gives one equation for each beam. 
As in the multislice derivation, we also neglect the second derivatives $d^2 \psi_{\bm{g}} / dz^2$, and obtain a set of coupled first-order differential equations for the beam amplitudes:

```{math}
:label: eq:bloch_coupled
\frac{d \psi_{\bm{g}}(z)}{dz} 
= 2 \pi i \left[ s_{\bm{g}} \psi_{\bm{g}}(z) + \sum_{\bm{h}} \frac{U_{\bm{g}-\bm{h}}}{2 k_0} \psi_{\bm{h}}(z) \right],
```
where 

```{math}
s_{\bm{g}} = \frac{k_0^2 - |\bm{k}_0 + \bm{g}|^2}{2 k_0}
```
is the excitation error of beam $\bm{g}$, which measures how far the reciprocal lattice point $\bm{g}$ lies from the Ewald sphere.
We can write these equations compactly in matrix form as

```{math}
:label: eq:bloch_ode
\frac{d \bm{\Psi}(z)}{dz} = 2 \pi i \, \bm{A} \, \bm{\Psi}(z),
```
where $\bm{\Psi}(z)$ is a vector containing all of the beam amplitudes $\psi_{\bm{g}}(z)$, and the structure matrix $\bm{A}$ has diagonal elements $A_{\bm{g}\bm{g}} = s_{\bm{g}} + U_{\bm{0}}/(2 k_0)$ and off-diagonal elements $A_{\bm{g}\bm{h}} = U_{\bm{g}-\bm{h}}/(2 k_0)$.

[](#eq:bloch_ode) is a first-order matrix differential equation, and this is the key to understanding the Bloch wave method. Its solution for a crystal of thickness $z$ is

```{math}
:label: eq:bloch_solution
\bm{\Psi}(z) = \exp\left(2 \pi i \bm{A} z \right) \bm{\Psi}(0).
```
To evaluate the matrix exponential, we compute the eigenvalue decomposition $\bm{A} = \bm{C} \bm{\Gamma} \bm{C}^{-1}$, where $\bm{\Gamma}$ is a diagonal matrix containing the eigenvalues $\gamma_j$, and the columns of $\bm{C}$ are the eigenvectors with elements $c_{\bm{g},j}$. 
The matrix exponential is then $\exp(2 \pi i \bm{A} z) = \bm{C} \exp(2 \pi i \bm{\Gamma} z) \bm{C}^{-1}$, where $\exp(2 \pi i \bm{\Gamma} z)$ is a diagonal matrix with elements $e^{2\pi i \gamma_j z}$. 
For a real crystal potential, $U_{\bm{g}-\bm{h}} = U^*_{\bm{h}-\bm{g}}$, the structure matrix is Hermitian, and the inverse $\bm{C}^{-1} = \bm{C}^\dagger$ is simply the conjugate transpose of $\bm{C}$. 
Each eigenvector is one Bloch wave: $c_{\bm{g},j}$ are the plane wave coefficients in [](#eq:bloch_wave), and the eigenvalue $\gamma_j$ gives the Bloch wavevector $\bm{k}_j = \bm{k}_0 + \gamma_j \hat{\bm{z}}$, which differs from the incident wavevector only along the beam direction.
Once $\bm{A}$ has been diagonalized, the wavefunction at any depth is given by the simple phase factors $e^{2\pi i \gamma_j z}$, as described in the next steps.

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
where $\psi_{\bm{g}}(0)$ represents the Fourier components of the input wavefunction. For plane wave illumination along $\bm{k}_0$, $\psi_{\bm{0}}(0) = 1$ and all other $\psi_{\bm{g}}(0) = 0$.

The expansion of $\psi_0(\bm{r})$ into Bloch waves is achieved by determining the weighting coefficients $\alpha_j$, which describe the contribution of each Bloch wave $b_j(\bm{r})$ to the input wavefunction. These coefficients are the elements of the vector $\bm{C}^{-1} \bm{\Psi}(0) = \bm{C}^{\dagger} \bm{\Psi}(0)$, or

```{math}
\alpha_j = \sum_{\bm{g}} c_{\bm{g},j}^* \psi_{\bm{g}}(0),
```
which relates the Bloch wave expansion coefficients $\alpha_j$ to the plane wave coefficients $\psi_{\bm{g}}(0)$ of the input wavefunction and the coupling coefficients $c_{\bm{g},j}$, which describe the relationship between the plane wave and Bloch wave bases.

#### 6. Output Wavefunction

At the exit surface ($z = z_{\text{max}}$, corresponding to the crystal thickness $t$), the real-space wavefunction $\psi(\bm{r})$ is expressed using the Bloch wave expansion. Dropping the common factor $e^{2\pi i \bm{k}_0 \cdot \bm{r}}$, as we did for the multislice method, the wavefunction is given by:

```{math}
\psi(\bm{r}) = \sum_{j} \alpha_j \left( \sum_{\bm{g}} c_{\bm{g},j} e^{2\pi i \bm{g} \cdot \bm{r}} \right) e^{2\pi i \gamma_j t}.
```
where $\alpha_j$ are the weighting coefficients determined by the input wavefunction, and $c_{\bm{g},j}$ are the Bloch wave coefficients that describe the periodic components of the wavefunction. The plane waves $e^{2\pi i \bm{g} \cdot \bm{r}}$ correspond to reciprocal lattice vectors $\bm{g}$, following the lattice periodicity.
The term $e^{2\pi i \gamma_j t}$ accounts for the propagation of each Bloch wave through the crystal thickness $t$, with $\gamma_j$ representing the eigenvalues that describe the wavevector components along $z$.

Note in particular how the thickness $t$ modulates the contribution of each Bloch wave to the final wavefunction at the exit surface, and how the Bloch wave method can be used to quickly compute output wavefunctions for multiple crystal thicknesses.
