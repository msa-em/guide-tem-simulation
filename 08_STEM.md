---
title: Scanning Transmission Electron Microscopy
numbering:
  enumerator: 1.%s
---

(stem_sims)=
## STEM Simulations with Converged Probes

(initial-conditions-of-the-electron-wave)=
### Initial Conditions of the Electron Wave

In STEM the objective lens acts on the electron beam before the beam reaches the specimen, hence the effects of the lens aberrations and aperture is incorporated in the initial conditions of the wavefunction.

The incident probe wavefunction is easiest to define in Fourier space

```{math}
:label: eq:fourier_probe

\begin{align} 
    \Psi_0(\bm{k}) = A(\bm{k}) \mathrm{e}^{-i\chi(\bm{k})},
\end{align}
```

where the amplitude, {math}`A(\bm{k})`, is the aperture function and the phase, {math}`\chi(\bm{k})`, is the aberration function introduced in section [%s](). The aperture is usually a disk with a radius given by the cutoff semiangle

```{math}
\begin{align}
    A(\bm{k}) = 
    \begin{cases}
    1, &  \lambda k = \alpha < \alpha_{cutoff} \\
    0,              & \text{otherwise} .
    \end{cases}
\end{align}
```

While the above definition is typical for most STEM experiments, the initial wavefunction can be defined using any other complex function and other initial conditions may be required for simulating more exotic experiments such phase plate STEM.

The probe is transferred to the specimen using an inverse Fourier transform and, in the same step, we can shift the by the Fourier shift theorem

```{math}
:label: eq:realspace_probe

\begin{align} 
    \psi_0(\bm{r}, \bm{r}_p) = \mathcal{F}_{\bm{q}}^{-1} \left[\mathrm{e^{-i\chi(\bm{k}) - 2 \pi i \bm{k}\cdot \bm{r}_p }} A(\bm{k}) \right] ,
\end{align}
```

where, {math}`\bm{r}_p`, is a specified probe position in relative to the atomic coordinates.

In Vis.[](#fig_probe) , we present an interactive visualization for exploring the relationship between size and shape of the incident probe and some parameters of the aperture and aberration functions.



```{figure} #app:probe_size
:name: fig_probe
**Left** The Fourier space initial wavefunction in equation [](#eq:fourier_probe). **Middle** The real space probe at the specimen in equation [](#eq:realspace_probe). **Right**  The real space intensity. Start by considering the effect of changing just the aperture and energy. In real space, the probe forms the diffraction-limited Airy disk pattern; increasing the aperture or energy decreases the radius of the pattern. Consider the differences between the STEM and SEM probes.  Next, add defocus. In Fourier space, the phase starts oscillating radially with a linearly decreasing period, the resulting probe grows and its radial intensity profile may have multiple peaks and valleys. Now try to decrease the aperture to include only the inner slowly varying part of the phase; the result is a smaller, more well-behaved probe. Add some spherical aberration and try to compensate by adding some defocus to flatten the phase inside the aperture resulting in a better probe (uncorrected STEM at Scherzer). Lastly, make a large probe and observe the diffraction fringes from self-interaction(boundary artifacts). This is the issue that should be avoided by increasing the size of the unit cell, as described in section [](). 
```

<!-- :::{figure} figures/probe.png
:name: vis:probe

[(Link to visualization)](https://boiling-wildwood-85903.herokuapp.com/voila/render/probes.ipynb) **Left** The Fourier space initial wavefunction in [](#eq:fourier_probe). **Right** The real space probe at the specimen in [](#eq:realspace_probe). Start by considering the effect of changing just the aperture and energy. In real space, the probe forms the diffraction-limited Airy disk pattern; increasing the aperture or energy decreases the radius of the pattern. Notice that changing the energy shifts the maximum of the axes in Fourier space. This is because the simulation uses a fixed wavefunction sampling; hence the maximum simulated scattering angle is decreased when the energy is increased according to Eq. [%s](). Next, add some defocus. In Fourier space, the phase starts oscillating radially with a linearly decreasing period, the resulting probe grows and its radial intensity profile may have multiple peaks and valleys. Now try to decrease the aperture to include only the inner slowly varying part of the phase; the result is a smaller, more well-behaved probe. Next, add some spherical aberration and try to compensate by adding some defocus to flatten the phase inside the aperture resulting in a better probe. Lastly, make a large probe and observe the diffraction fringes from self-interaction. This is the issue that should be avoided by increasing the size of the unit cell, as described in section [](wraparound).
::: -->

:::{figure} figures/probe_overlap.png
:name: fig_probe_overlap
:align: center

**Estimating probe wraparound errors.** - [link to notebook](https://github.com/tem-elements/tem-elements/blob/main/notebooks/Probe_overlap.ipynb)
:::

Because STEM simulations can require long computation times, we often try to reduce the size of the simulation cell as much as possible. However, we must be careful to use a simulation cell large enough to hold the STEM probe. In particular, we need to consider the size of the probe throughout the full simulation cell volume. For an empty cell, the probe will have a maximum size at either the entrance or exit surface, depending on the probe defocus.

[](#fig_probe_overlap) shows an interactive demo for testing the overlap of a STEM probe for different microscope parameters and cell dimensions. Try setting ... Note that this example only considers an empty cell volume. Atoms inside the simulation volume will scatter the electron beam, with heavier elements scattering electrons to higher angles. We recommend positioning individual STEM probes at the positions where the highest scattering is expected (for example directly on or adjacent to the thickest atomic columns), and carefully checking the dimensions of the probe at the exit surface. This is especially important for PRISM simulations, as the cropping box around the STEM probe can be significantly smaller than the full cell dimensions for high PRISM interpolation factors {cite:t}`ophus_fast_2017`.

(imaging)=
### Imaging

The initial wavefunction is passed through the specimen potential according to the multislice algorithm. The exit wavefunction, {math}`\psi_t(\bm{r}, \bm{r}_p)`, is then diffracted to the detector plane using another Fourier transform

```{math}
\begin{align*}
    \Psi_t(\bm{k}, \bm{r}_p) = \mathcal{F}_{\bm{r}} [\psi_t(\bm{r}, \bm{r}_p)] .
\end{align*}
```

The square modulus, {math}`|\Psi_t(\bm{k}, \bm{r}_p)|^2` is the diffraction pattern for the probe position {math}`\bm{r}_p`.

Every STEM mode requires calculating a diffraction pattern for each shifted initial wavefunction in the scan region, the modes differ only by how the detector geometry is applied to the diffraction patterns.

In bright-field and annular-dark-field STEM the detector integrates the CBED pattern on a region in the diffraction plane. This is equivalent to multiplying with a detector function

```{math}
\begin{align*}
    g(\bm{r}_p) = \int\int D(\bm{k})|\Psi_t(\bm{k}, \bm{r}_p)|^2 ~d^2 \bm{k}
\end{align*}
```


where

```{math}
\begin{align*}
    D(\bm{k}) = 
    \begin{cases}
    1, &  \alpha_a < \lambda k < \alpha_b \\
    0,              & \text{otherwise}
    \end{cases}
\end{align*}
```

If {math}`D(\bm{k})` is a small point on the axis then the measurement is a bright field image. If {math}`D(\bm{k})` is a large annulus covering the high angle scattering then the measurement is an annular dark field image.

The image, {math}`g(\bm{r}_p)`, is the collection of integrated intensities at all sample positions in a (typically) rectangular region. The scan region may be chosen independently of the supercell and may cover only part of the super cell. In the case of a super cell consisting smaller periodic units, computation can be saved by choosing a scan window covering just one of the periodic units.

:::{figure} figures/annular_integrals.png
:name: vis:annular
:alt: (Link to visualization)

[(Link to visualization)](https://boiling-wildwood-85903.herokuapp.com/voila/render/annular_detector.ipynb)
:::

%  We can significantly limit the computational cost, by choosing the largest possible spacing between the sample positions, or probe step size, $\Delta r_p$. 

%  \begin{align*}

%      \Delta r_p = 0.9 \frac{1}{4 \lambda \alpha_{cutoff}}.

%  \end{align*}

%  It is important to remember the difference between the wavefunction sampling and probe step, both in units of Angstrom. The probe step only refers to the spacing of the initial probe, and it is entirely independent from the grid used to sample the wavefunctions and potentials.

%  In Vis. \ref{}, we present a visualization for exploring how the integration region of the detector influence the image.

(differential-phase-contrast)=
### Differential Phase Contrast
The phase problem, namely the loss of phase information when taking a measurement, is well known in many fields including electron microscopy.  Because detectors collect the square modulus of the exit wave ({math}`|\Psi_t(\bm{k}, \bm{r}_p)|^2`), much of the phase information is lost. Methods that can reconstruct the phase of the sample provide dose-efficient imaging of materials and the ability to simultaneously characterize heavy and light elements. There are a number of approaches to recovering the phase of the sample, some of which will be discussed in the [4D-STEM](#id-4d-stem) section.

When an electron probe interacts with a sample potential the center of mass (CoM) of the beam changes. Depending on the sample potential and the size of the probe, this will take the form of either a rigid disk shift, for low-frequency features, or a change in the distribution of signal within a disk, for high-frequency features, such as atoms {cite:p}`cao2018theory`. 

First proposed by {cite:t}`dekkers1974differential` replacing conventional detectors with segemented detectors can caputre this change in center of mass. Differential phase contrast (DPC) measurements are made by differentiating signal from oppposite segments, and this signal is proportional to the phase of the sample. Fig. illustrates this DPC approach with the mixed nanoparticle on carbon sample, and compared to the bright-field and dark-field images, the signal from both the heavy and light components is more clear. The contrast transfer function for DPC peaks for the convergence angle of the probe, so simulations with these detectors can be helpful for optimizing experimental parameters.

(id-4d-stem)=
### 4D-STEM (CO)


