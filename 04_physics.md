---
title: Physics Concepts
numbering:
  enumerator: 1.%s
math:
  \ii: '{i\mkern1mu}'
  \invFFT: \mathcal{F}^{-1}_{\mathbf{k}\to\mathbf{r}}
  \FFT: \mathcal{F}_{\mathbf{r}\to\mathbf{k}}
  \angstroms: \text{\normalfont\AA}
bibliography:   
  - refs
github: https://github.com/msa-eom/guide-tem-simulation
---

(physical-concepts-ts)=
# Physical Concepts (TS)

The properties of atoms, molecules and solids are fundamentally determined by quantum mechanics. In this modern theory of physics, some classical concepts such as the electrostatic potential and the propagation of waves are carried over essentially unchanged, whereas many familiar intuitions fail on the level of the very small. In particular, quantum objects including electrons exhibit both wave and particle properties depending on how they are observed — in the context of electron scattering, matter-wave interference is of central importance. Mathematically, both free propagating electrons and those bound into atoms are described mathematically as complex waves, via so-called electron wavefunctions.

(electron-wavefunctions-ts)=
## Electron Wavefunctions

An electron wavefunction, typically denoted by $\psi$, is a mathematical description of the quantum state of an isolated quantum system. The wavefunctions are complex-valued, and thus not directly observable, but their squared amplitudes can be interpreted as probabilities to find the electrons in particular states that do correspond to physically observable properties of the system. Wavefunctions naturally live in a many-dimensional mathematical space called the Hilbert space. By choosing a basis of representation, they can be represented in real space, typically using either Cartesian (convenient for planewave propagation) or spherical coordinates (convenient for scattering).

Two kinds of wavefunctions are commonly encountered in the context of electron scattering: freely propagating plane and spherical waves, and wavefunctions of electrons bound by (a collection of) atoms. Plane-waves can be mathematically denoted as

```{math}
:label: eq:plane-waves
\psi(\bm{r}, t)=\frac{1}{\sqrt{V}} \mathrm{e}^{i\left(\bm{k} \cdot \bm{r}-\omega t\right)},
```

where $\bm{r}$ is the three-dimensional position vector, $t$ is the time, $\bm{k}$ with magnitude $\left|\bm{k}\right|= k = 2 \pi/\lambda$ is the three-dimensional wavevector for wavelength $\lambda$, and $\omega$ is the angular frequency of the wave. $V$ is a volume factor that ensure that integral of the wavefunction over all of space is equal to unity.

Spherical waves can be described similarly,

```{math}
:label: eq:spherical-waves
\psi(\bm{r}, t)=\frac{1}{\sqrt{V}} \frac{\mathrm{e}^{i\left(\bm{k}\left|\bm{r}-\bm{r}^{\prime}\right|-\omega t\right)}}{\left|\bm{r}-\bm{r}^{\prime}\right|},
```

where $\bm{r}^{\prime}$ now denotes the location of a scattering center, and the denominator $\left|\bm{r}-\bm{r}^{\prime}\right|$ ensures the correct decay of the amplitude as a function of inverse distance.

To describe the propagation of these waves in free space, the imaginary term in the exponential describes periodic variation. Thus for any fixed sum of the spatial and temporal terms in the exponent, the wavefunction has the same value; these represent the planes (or shells) of constant amplitude. To calculate the wave in an arbitrary position, we can simply substitute the new spatial location and time arguments to calculate the resulting amplitude.

(bound-systems-ts)=
## Bound Systems: The Schrödinger Equation (TS)

To derive the wavefunction of a bound system, we need to solve the quantum "equation of motion" for the electron(s) in the corresponding confining potential: this is the Schrödinger equation. The Schrödinger equation gives the fundamental mathematical description of quantum systems, describing the time-evolution of their wavefunction. The equation for a single non-relativistic particle can in the position representation be written as

```{math}
:label: eq:Schrodinger_time
    \hbar \frac{\partial}{\partial t} \psi(\bm{r}, t)=\left[-\frac{\hbar^{2}}{2 m} \frac{\partial^{2}}{\partial \bm{r}^{2}}+V(\bm{r}, t)\right] \psi(\bm{r}, t),
```

where $\hbar$ is the reduced Planck constant, $m$ is the electron mass, and $V(\bm{r}, t)$ is the potential, e.g., Coulomb potential of a nucleus in the case of an atom. Mathematically, the Schrödinger equation is linear partial differential equation that assigns a complex number to each point $\bm{r}$ at time $t$, which can be interpreted as the probability amplitude of the electron.

The equation can only be solved analytically for a handful of simple cases, of which the hydrogen atom is an illustrative example. For many such cases, we look for stationary solutions and can thus use the time-independent form the Schrödinger equation, written as

```{math}
:label: eq:Schrodinger
    \left[-\frac{\hbar^{2}}{2 m} \frac{\partial^{2}}{\partial \bm{r}^{2}}+V(\bm{r})\right] \psi(\bm{r}) = E\psi(\bm{r}),
```

where wavefunctions $\psi(\bm{r})$ are the eigenvectors and energies $E$ the eigenvalues of the system.

(electrostatic-potentials-ts)=
## Electrostatic Potentials (TS)

The electrostatic potential of a specimen determines not only how the electrons of the system are bound, but also how transmitting electrons scatter via the Lorentz force. It therefore connects the properties of the material to the resulting images or diffraction patterns. The electrostatic potential is fundamentally speaking derived from the electron density of the atoms in a specimen, which is described by their quantum mechanical many-body wavefunction. However, this is only rarely analytically solvable, and various approximations may be needed.

In the case of the hydrogen atom, an analytical solution is possible. To describe the confining effect of the proton of the nucleus on the lone valence electron, we need to include the attractive Coulomb potential in the Hamiltonian. Using the time-independent Schrödinger equation, this can be written in spherical coordinates as

```{math}
:label: eq:H_Schrodinger
\left(-\frac{\hbar^{2}}{2 m} \frac{\partial^{2}}{\partial \bm{r}^{2}}-\frac{q_e^2}{4 \pi \varepsilon_{0} r}\right) \psi(r, \theta, \varphi)=E \psi(r, \theta, \varphi),
```

where $q_e$ is the elementary charge, $\varepsilon_0$ is the permittivity of free space, $r$ the distance from the nucleus, and $\theta$ and $\varphi$ are respectively the spherical polar and azimuthal angles. The well-known normalized solution can be expressed as

```{math}
:label: eq:H_solution
\psi_{n \ell m}(r, \theta, \varphi)=\sqrt{\left(\frac{2}{n a_0}\right)^{3} \frac{(n-\ell-1) !}{2 n(n+\ell) !}} \mathrm{e^{-\rho / 2}} \rho^{\ell} L_{n-\ell-1}^{2 \ell+1}(\rho) Y_{\ell}^{m}(\theta, \varphi),
```

where the quantum numbers ($n$ $\ell$ $m$) numerate the possible states of the electron, $\rho=\frac{2 r}{n a_0}$ with the Bohr radius $a_0=\frac{4 \pi \varepsilon_{0} \hbar^{2}}{m q_e^2}$, and $L_{n-\ell-1}^{2 \ell+1}(\rho)$ is a generalized Laguerre polynomial of degree $n-\ell-1$, and $Y_{\ell}^{m}(\theta, \varphi)$ is a spherical harmonic function of degree $\ell$ and order $m$.

To calculate the electron density, we need to select a set of quantum numbers ($n$ $\ell$ $m$); for the ground state, we can choose $n = 1$, $\ell = m = 0$. The wavefunction thus simplifies to

```{math}
:label: eq:H_wavefunction
\psi_{1 0 0}(r) = \frac{1}{\sqrt{\pi} a_0^{3 / 2}} \mathrm{e^{-r / a_{0}}},
```

whose squared norm gives the (non-relativistic) electron density of hydrogen as

```{math}
:label: eq:H_density
\rho(r) = |\psi_{1 0 0}(r)|^2 = \frac{1}{\pi a_0^3} \mathrm{e^{-2 r / a_0}}.
```

(scattering-from-a-potential-ts)=
## Scattering from a Potential

To understand how electrons scatter from a potential, we can consider an electron plane wave incident on an isolated atom. The wave interacts with the electrostatic potential of the nucleus and electrons of the atom, and an outgoing spherical wave is generated. To understand the effect of the atom on the wave, we need to calculated the distribution of scattered intensity, which is not isotropic due to the initial linear momentum of the incident wave.

The scattering problem can be formally treated by finding a solution of the Schrödinger equation [](#eq:Schrodinger) for the incident electron inside the scattering atom (with electron coordinates $\bm{r^{\prime}}$)

```{math}
    \left[-\frac{\hbar^{2}}{2 m} \frac{\partial^{2}}{\partial \bm{r^{\prime \textmd{2}}}}+V(\bm{r^{\prime}})\right] \psi(\bm{r^{\prime}}) = E\psi(\bm{r^{\prime}}),
```

which we can re-write as

```{math}
    \left(\nabla^{2}+k_{0}^{2}\right) \psi\left(\bm{r}^{\prime}\right)=U\left(\bm{r}^{\prime}\right) \psi\left(\bm{r}^{\prime}\right)
```

with the substitutions 

```{math}
\begin{aligned}
    k_{0}^{2} & \equiv \frac{2 m E}{\hbar^{2}}, \\ U\left(\bm{r}^{\prime}\right) & \equiv \frac{2 m V\left(\bm{r}^{\prime}\right)}{\hbar^{2}}.
\end{aligned}
```

This can be formally solved with the help of Green's function $G(\bm{r},\bm{r^{\prime}})$ (we refer the reader to {cite:t}`fultz_transmission_2013` for detail) to yield a sum of the incident and scattered components

```{math}
:label: eq:scattering_schrodinger
    \psi(\bm{r}) = \psi_{\mathrm{inc}}(\bm{r}) + \psi_{\mathrm{scatt}}(\bm{r}) =
    \mathrm{e}^{\mathrm{i} k_{0} \cdot \bm{r}}+\frac{2 m}{\hbar^{2}} \int V\left(\bm{r}^{\prime}\right) \psi\left(\bm{r}^{\prime}\right) G\left(\bm{r}, \bm{r}^{\prime}\right) \mathrm{d}^{3} \bm{r}^{\prime}.
```

While this formal solution is in principle exact, it is in the form of an implicit integral equation whose solution $\psi$ appears both inside and outside the integral – indeed, we have merely transformed the differential Schrödinger equation into an integral equation without solving anything. To proceed, some approximation is needed.

(born-approximation-ts)=
### Born approximation for electron scattering

A widely used approximation to solve [](#eq:scattering_schrodinger) is the first Born approximation, whereby we replace the full $\psi$ within the integral simply by the incident planewave. In effect, this approximation makes the assumption that the incident wave is not diminished and scattered only once by the material, which is valid when scattering is weak.

By further assuming that the detector is far from the scatterer, which allows us to work with outgoing planewaves instead of spherical waves, and by aligning the outgoing wavevector $\bm{k}$ with the direction $(\bm{r} - \bm{r^{\prime}})$, as well as assuming that the origin of the coordinate system is near the atom so that $|\bm{r}| \gg |\bm{r^{\prime}}|$, we can write the Green's function as

```{math}
    G\left(\bm{r}, \bm{r}^{\prime}\right) \simeq-\frac{1}{4 \pi} \frac{\mathrm{e}^{i \bm{k} \cdot\left(\bm{r}-\bm{r}^{\prime}\right)}}{|\bm{r}|}.
```

By substituting this to [](#eq:scattering_schrodinger) we get

```{math}
\begin{align}
\psi(\bm{r}) &\simeq \mathrm{e}^{\mathrm{i} k_{0} \cdot \bm{r}}-\frac{m}{2 \pi \hbar^{2}} \int V\left(\bm{r}^{\prime}\right) \mathrm{e}^{\mathrm{i} k_{0} \cdot \bm{r}^{\prime}} \frac{\mathrm{e}^{\mathrm{i} k \cdot\left(\bm{r}-\bm{r}^{\prime}\right)}}{|\bm{r}|} \mathrm{d}^{3} \bm{r}^{\prime}\\
&=\mathrm{e}^{\mathrm{i} k_{0} \cdot \bm{r}}-\frac{m}{2 \pi \hbar^{2}} \frac{\mathrm{e}^{\mathrm{i} \bm{k} \cdot \bm{r}}}{|\bm{r}|} \int V\left(\bm{r}^{\prime}\right) \mathrm{e}^{\mathrm{i}\left(\bm{k}_{0}-\bm{k}\right) \cdot \bm{r}^{\prime}} \mathrm{d}^{3} \bm{r}^{\prime}.
\end{align}
```

By further defining $\Delta \bm{k} \equiv \bm{k}-\bm{k}_{0}$, we can write the scattered part of the wave as

```{math}
\psi_{\mathrm{scatt}}(\Delta \bm{k}, \bm{r})=\frac{\mathrm{e}^{\mathrm{i} k \cdot \bm{r}}}{|\bm{r}|} f(\Delta \bm{k}),
```

which is a function of only the difference between the incident and outgoing wavevectors $\Delta \bm{k}$, describes the change in the direction of the scattered radiation with respect to the incident direction (i.e. momentumn transfer). The factor

```{math}
:label: eq:formfactor
f(\Delta \bm{k}) \equiv-\frac{m}{2 \pi \hbar^{2}} \int V\left(\bm{r}^{\prime}\right) \mathrm{e}^{-\mathrm{i} \Delta \bm{k} \cdot \bm{r}^{\prime}} \mathrm{d}^{3} \bm{r}^{\prime},
```

is called the atomic form factor (or electron scattering factor {cite:p}`kirkland_advanced_2010` and it describes the angular distribution of scattered intensity. In the first Born approximation it corresponds to the Fourier transform of the scattering potential ($f(\Delta \bm{k})=\mathcal{F}_k[V(r)]$).

(atomic-form-factors-ts)=
### Atomic form factors (TS)

The atomic form factors thus describe the angular amplitude for scattering of a single electron of by a single atom. While the first Born approximation is inadequate for describing real specimens (as electrons typically scatter multiple times when passing through a crystal), this description is quite useful since it relates the three-dimensional Fourier transform of the atomic potential to the scattering amplitude.

In the general case, we can describe the potential as a sum of a negative term due to the Coulomb attraction of the nucleus with atomic number $Z$ and a positive term arising from the Coulomb repulsion of the atomic electrons with electron density $\rho(\bm{r})$

```{math}
V(\bm{r})=-\frac{Z e^{2}}{|\bm{r}|}+\int_{-\infty}^{+\infty} \frac{e^{2} \rho\left(\bm{r}^{\prime}\right)}{\left|\bm{r}-\bm{r}^{\prime}\right|} \mathrm{d}^{3} \bm{r}^{\prime}.
```

By substituting this into [](#eq:formfactor) and defining a new variable $\bm{R} \equiv \bm{r} - \bm{r}^{\prime}$, so that $\bm{r} = \bm{R} - \bm{r}^{\prime}$ and rearranging, we get

```{math}
f(\Delta \bm{k})=\frac{m Z e^{2}}{2 \pi \hbar^{2}} \int_{-\infty}^{+\infty} \frac{1}{|\bm{r}|} \mathrm{e}^{-\mathrm{i} \Delta k \cdot \bm{r}} \mathrm{d}^{3} \bm{r} -\frac{m e^{2}}{2 \pi \hbar^{2}} \int_{-\infty}^{+\infty} \frac{1}{|\bm{R}|} \mathrm{e}^{-\mathrm{i} \Delta \bm{k} \cdot \bm{R}} \mathrm{d}^{3} \bm{R} \int_{-\infty}^{+\infty} \rho\left(\bm{r}^{\prime}\right) \mathrm{e}^{-\mathrm{i} \Delta \bm{k} \cdot \bm{r}^{\prime}} \mathrm{d}^{3} \bm{r}^{\prime}.
```

The two first integrals are simply Fourier transforms of $1/r$ that each yield $4 \pi / \Delta k^2$, giving the general expression for the electron scattering factor of an atom as

```{math}
f(\Delta \bm{k})=\frac{2 m e^{2}}{\hbar^{2} \Delta k^{2}}\left(Z-\int_{-\infty}^{+\infty} \rho(\bm{r}) \mathrm{e}^{-\mathrm{i} \Delta \bm{k} \cdot \bm{r}} \mathrm{d}^{3} \bm{r}\right).
```

For kinematic scattering (where diffraction intensity $I(k)=|\psi(k)|^{2}\approx\left|f(k)\right|^{2}$) and within the first Born approximation for electrons (i.e. electron scattering factor is the Fourier transform of the total electrostatic potential) this gives the Mott-Bethe formula for a single atom,

```{math}
f(\Delta \bm{k})=\frac{2 m e^2}{\hbar^{2} \Delta k^{2}} \left(Z-\frac{m c^2}{e^2}f_{x}(\Delta \bm{k})\right),
```

relating electron scattering factors with x-ray scattering factors $f_{x}(k)=\mathcal{F}_k[\rho(r)]$, where $\rho(r)$ is the electron density. This comparison highlights that unlike x-ray scattering (where the nucleus is too heavy to accelerate in an interaction with a nearly massless photon), electron scattering is sensitive to both nuclear and electron charges, with the latter contribution being relatively enhanced for small electron scattering angles (corresponding to small momentum vectors $k$). The Mott-Bethe formula is a convenient way to obtain the potential distribution from a charge distribution including the nuclear point charge, effectively by solving Poisson's equation in reciprocal space. 

As an example, we can calculate the x-ray form factor for hydrogen in spherical coordinates (c.f. [](#eq:H_density), assuming that the incident electron wavevector $\bm{k} = (0, 0, k)$ is parallel to the $z$ axis. In this spherically symmetric case, the x-ray form factor reduces to a radial integral

```{math}
f_\mathrm{x}(\Delta k)=\frac{e^2}{m c^2} \int_{0}^{\infty} \rho(r) \sin (2 \pi \Delta k\, r)\, r\, \mathrm{d}r,
```

where the exponential term of the Fourier transform has been expanded as a trigonometric integral. We can calculate the integral explicitly as

```{math}
f_{x}(\Delta k) = \frac{e^2}{m c^2} \int_{0}^{\infty} \frac{r \mathrm{e}^{-2r / a_0}}{\pi a_0^{3}} \sin (2 \pi \Delta k \,r)\, \mathrm{d} r \nonumber
     = \frac{e^2}{m c^2}\frac{1}{\left(1+\pi^2 \Delta k^2 a_0^{2}\right)^{2}}.
```

Using the definition of the Bohr radius ($a_{0}=4 \pi \varepsilon_{0} \hbar^{2} / m e^{2}$), the electron scattering factor for hydrogen ($Z = 1$) can then further be written via by the Mott-Bethe formula as 

```{math}
f(\Delta k)=\frac{a_0}{2 \varepsilon_0 \Delta k^{2}} \left(1-\frac{1}{\left(1+\pi^2 \Delta k^2 a_0^{2}\right)^{2}}\right).
```

Hydrogen is the only case that is analytically solvable; in general neither the density nor the atomic form factor can be directly written down. Instead, these have been calculated using various approximations to the true multielectron wavefunctions, and tabulated for isolated atoms of all elements as isolated atomic potentials.

(isolated-atomic-potentials-ts)=
## Isolated Atomic Potentials (TS)

Since the Schrödinger equation cannot be solved analytically even for most molecules, let alone solid-state systems, several kinds of approaches have been developed to obtain approximate solutions. These accurate but computationally very expensive techniques have been used to parametrize what are called isolated atomic potentials – or, equivalently in reciprocal space, electron scattering factors – which describe the potential of a specimen as a sum of isolated, non-interacting atom potentials. This approximation is often called the independent atom model.

A potential parametrization is a numerical fit to such first principles calculations of electron atomic form factors that describe the radial dependence of the potential for each element. One of the most widely used parametrizations is the one published in {cite:t}`kirkland_advanced_2010`, who fitted Dirac-Fock scattering factors with combination of Gaussians and Lorentzians. In 2014, Lobato and Van Dyck improved the quality of the fit further, using hydrogen's analytical non-relativistic electron scattering factors as basis functions to enable the correct inclusion of all physical constraints {cite:p}`lobato_accurate_2014`.

An example of independent atom model scattering factors and potentials for several elements up to $Z = 32$ is shown in [](#fig:atomic_potentials).

```{figure} #parametrized_potentials
:name: fig:atomic_potentials
Independent atom model potentials and scattering factors for several elements.
```

(dft-potentials-ts)=
## Density Functional Theory Potentials (TS)

Since the complicated many-body interactions of multiple electrons mean that wavefunction cannot in general be analytically solved, further approximations are needed. Density functional theory (DFT) is the most prominent one, and is widely use for modeling the electronic structure of molecules and solids. In DFT, the combinatorially intractable many-body problem of $N$ electrons with $3N$ spatial coordinates is reduced to a solution for the three spatial coordinates of the electron density that can be variationally reached. This approximation would in principle be exact, but a term that describes electron exchange and correlation is not analytically known and must be approximated.

Although the ground-state electron density for all electrons, including those in the core levels and in the valence, can in principle be solved within the DFT framework, the electron wavefunctions rapidly oscillate near the nuclei, making a numerical description computationally very expensive. To make calculations practical, some partition of the treatment of the cores and the valence is therefore typically needed. Pseudopotential {cite:p}`schwerdtfeger_pseudopotential_2011` and projector-augmented wave (PAW) methods {cite:p}`blochl_projector_1994` have in recent years matured to offer excellent computational efficiency. The core electrons are not described explicitly in either method, but in the former are replaced by a smooth pseudo-density, while in the latter, by smooth analytical projector functions in the core region.

Inverting the projector functions allows the exact core electron density to be analytically calculated, making the PAW method arguably ideally suited for efficient and accurate *ab initio* all-electron electrostatic potentials. The PAW method is accordingly the approach chosen for *ab*TEM, specifically via the grid-based DFT code GPAW (more details on the method can be found in the literature {cite:p}`blochl_projector_1994,enkovaara_electronic_2010`. Notably for our purposes here, GPAW allows the full electrostatic potentials (with smeared nuclear contributions) to be efficiently calculated for molecules and solids, and seamlessly used for electron scattering simulations using *ab*TEM.

To come back to our hydrogen example, we can analytically calculate its full electrostatic potential by solving the Poisson equation for the charge density of the electron (charge $-q_e$) combined with the Coulomb potential of the nucleus ($+q_e$). The Poisson equation is

```{math}
\nabla^{2} V(r) = -\frac{4\pi}{\varepsilon_0} \left(-q_e \rho(r) + q_e \delta(r)\right),
```

with the solution

```{math}
V(r) = \frac{q_e}{4 \pi \epsilon_{0}} \frac{\mathrm{e^{-2 r / a_0}}}{r}\left(1 + \frac{r}{a_0}\right).
\label{eq:H_potential}
```

This expression is singular at the origin due to the point charge of the nucleus, but that singularity gets smeared out due and numerically discretized in practical calculations (by default, GPAW smears the nuclear potentials with a Gaussian function).

In [](#fig:H_atom), we plot the exact solution of [](#eq:H_potential) against the Kirkland IAM parameterization as well as a DFT calculation using the GPAW code. While all models agree perfectly near the nucleus (the discontinuous appearance of the line is due to a finite computational grid with a spacing of 0.01 Å), where the potential is the strongest, small differences emerge further away. In the case of the DFT model, the simulation box is finite in size, and thus the continuity requirement of the wavefunction and corresponding density slightly affect the long-range part of the calculated potential. The choice of exchange-correlation functional may also slightly influence the results.

```{figure} #H_potential_comparisons
:name: fig:H_atom
Comparing the radial electrostatic potential of the hydrogen atom exactly solved from the Schrödinger equation, parametrized by the independent atom model, and calculated by DFT.
```

Although IAM potentials are useful for many purposes, they do neglect chemical bonding, which may be measurable and of interest. To illustrate this difference, an interactive comparison of the IAM and DFT scattering potentials of the H$_2$ molecule at different distances between the H atoms is shown below.

```{figure} #H2_potential
:name: fig:H_atom
The difference between the indepdendent atom model potential to the DFT potential for the hydrogen molecule as a function of the distance between the H atoms.
```

(numerical-solutions-of-the-schrodinger-equation-co)=
### Numerical Solutions of the Schrödinger Equation (CO)

As discussed previously, the Schrödinger equation typically cannot be solved analytically in complex systems. Therefore, in order to perform electron scattering simulations, we must calculate numerical solutions of [](#eq:Schrodinger_time) for electron waves. First, we define the {cite:t}`debroglie1925recherches` wavelength of a free electrons (corrected for relativistic effects) as

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

where {math}`{\nabla_{xy}}^2 = \partial^2/\partial x^2 + \partial^2/\partial y^2`. This equation shows the overall numerical recipe we will use; when the wavefunction {math}`\psi_0(\bm{r})` is at position {math}`z_0`, we will evaluate the operators on the right hand side over a distance {math}`\Delta z` to calculate the new wavefunction {math}`\psi(\bm{r})` at position {math}`z_0 + \Delta z`. {cite:t}`kirkland2020` gives the formal operator solution to [](#eq:Shrodinger_electron) as

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

is a thin slice of the potential. Unfortunately, even with the above approximations, [](#eq:Shrodinger_simple) cannot be solved in closed form due to the two non-commuting operators. Instead, we will solve it numerically by using a split-step method, where we calculate solutions for each operator independently, and alternating their application to the electron wavefunction. This solution was introduced by {cite:t}`cowley1957scattering` and is known as the multislice method.
