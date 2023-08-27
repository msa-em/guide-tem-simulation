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

In [](#fig_probe) , we present an interactive visualization for exploring the relationship between size and shape of the incident probe and some parameters of the aperture and aberration functions.



```{figure} #app:probe_size
:name: fig_probe
**Left** The Fourier space initial wavefunction in equation [](#eq:fourier_probe). **Middle** The real space probe at the specimen in equation [](#eq:realspace_probe). **Right**  The real space intensity. Start by considering the effect of changing just the aperture and energy. In real space, the probe forms the diffraction-limited Airy disk pattern; increasing the aperture or energy decreases the radius of the pattern. Consider the differences between the STEM and SEM probes.  Next, add defocus. In Fourier space, the phase starts oscillating radially with a linearly decreasing period, the resulting probe grows and its radial intensity profile may have multiple peaks and valleys. Now try to decrease the aperture to include only the inner slowly varying part of the phase; the result is a smaller, more well-behaved probe. Add some spherical aberration and try to compensate by adding some defocus to flatten the phase inside the aperture resulting in a better probe (uncorrected STEM at Scherzer). Lastly, make a large probe and observe the diffraction fringes from self-interaction(boundary artifacts). This is the issue that should be avoided by increasing the size of the unit cell, as described in section [](). 
```


Because STEM simulations can require long computation times, we often try to reduce the size of the simulation cell as much as possible. However, we must be careful to use a simulation cell large enough to hold the STEM probe. In particular, we need to consider the size of the probe throughout the full simulation cell volume. For an empty cell, the probe will have a maximum size at either the entrance or exit surface, depending on the probe defocus.

Use the drop down menu in [](#fig_probe), to see what the probe looks like when it has wrap-around artifacts. Try modifying the aberrations and notice as the probe gets bigger the error gets worse. Note that this example only considers an empty cell volume. Atoms inside the simulation volume will scatter the electron beam, with heavier elements scattering electrons to higher angles. We recommend positioning individual STEM probes at the positions where the highest scattering is expected (for example directly on or adjacent to the thickest atomic columns), and carefully checking the dimensions of the probe at the exit surface. This is especially important for PRISM simulations, as the cropping box around the STEM probe can be significantly smaller than the full cell dimensions for high PRISM interpolation factors {cite:t}`ophus_fast_2017`.

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

We can significantly limit the computational cost, by choosing the largest possible spacing between the sample positions, or probe step size, $\Delta r_p$. 

\begin{align*}

      \Delta r_p = 0.9 \frac{1}{4 \lambda \alpha_{cutoff}}.

\end{align*}

It is important to remember the difference between the wavefunction sampling and probe step, both in units of Angstrom. The probe step only refers to the spacing of the initial probe, and it is entirely independent from the grid used to sample the wavefunctions and potentials.

In [](#fig_stem_processing), we present a visualization for exploring how the integration region of flexible circular and annular detectors influence the image contrast. More discussion of flexible detectors in STEM experiments is in [4D-STEM](#id-4d-stem). In bright field imaging, we generally use a detector slightly bigger than the probe size to integrate all the information from the directly transmitted beam. In annular dark field (ADF) imaging the contrast is approximately $Z^\alpha$, where Z is the atomic number and $alpha$ is anywhere between 1.3 and 2 depending on the detector geometry {cite:p}`treacy2011z`. In a high angle annular dark field (HAADF) experiment, the contrast is even more strongly dominated by the heavy atoms.  ADF imaging is a linear technqiue and relatively robust to sample thickness, but as can be oserved in [](#fig_stem_processing) is not well suited for visualizing light elements. Annular bright field imaging (ABF) uses a colleciton angle matching the outer ring of the bright field disk, and capture signals from both heavy and light elements {cite:p}`okunishi2009visualization`. 


```{figure} #app:stem_post_processing
:name: fig_stem_processing
Aluminum, iron and gold nanoparticles on a carbon film: **Left** image from a circular detector. **Middle** image from an annular detector **Right**  differential phase contrast reconstruction. Notice how changing the collection angles impacts the contrast and signal to noise. 
```

(differential-phase-contrast)=
### Differential Phase Contrast
The phase problem, namely the loss of phase information when taking a measurement, is well known in many fields including electron microscopy.  Because detectors collect the square modulus of the exit wave ({math}`|\Psi_t(\bm{k}, \bm{r}_p)|^2`), much of the phase information is lost. Methods that can reconstruct the phase of the sample provide dose-efficient imaging of materials and the ability to simultaneously characterize heavy and light elements. There are a number of approaches to recovering the phase of the sample, some of which will be discussed in the [4D-STEM](#id-4d-stem) section.

When an electron probe interacts with a sample potential the center of mass (CoM) of the beam changes. Depending on the sample potential and the size of the probe, this will take the form of either a rigid disk shift, for low-frequency features, or a change in the distribution of signal within a disk, for high-frequency features, such as atoms {cite:p}`cao2018theory`. 

First proposed by {cite:t}`dekkers1974differential` replacing conventional detectors with segemented detectors can capture this change in center of mass. Differential phase contrast (DPC) measurements are made by differentiating signal from oppposite segments, and this signal is proportional to the phase of the sample. [](#fig_stem_processing) illustrates a DPC approach for the mixed nanoparticle on carbon structure, and compared to the bright-field and dark-field images, the signal from both the heavy and light components is more clear. The DPC approach is also more efficient than the ABF mode. The contrast transfer function for DPC peaks for the convergence angle of the probe, so simulations with these detectors can be helpful for optimizing experimental parameters.

(id-4d-stem)=
### 4D-STEM 

In a conventional imaging experiment, detectors of fixed collection angle integrate the signal from the exit wave. By contrast, in a 4D-STEM experiment, at each probe position a full diffraction pattern is recorded. A 4D-dataset contains spatially resolved structural information, and post-processing can be used to recover features not possible with a conventional imaging technqiues. 

There are myriad possibilities with 4D-STEM, and we refer you to other references for a more complete discussion {cite:p}`ophus2019four`. It is difficult to visualize a 4D-dataset, so most reconstruction approaches map this 4D-dataset into real or reciprocal space reflections of the full dataset. One of the more straight forward and common applications of 4D-STEM virtual imaging. In a virtual imaging reconstruction, a detector is applied in reciprocal space to create a real space image, and we refer to the detector as virtual because it is applied computationally after the data is acquired. Virtual detectors can take any form including binary filters that match the geometry of conventional detectors and circular detectors to pick-up specific Bragg disks. Similiarly, a virtual mask can be applied in real space to show diffraction patterns from a subset of the field of interest.

```{figure} #app:stem_graphene
:name: fig_stem_graphene
Using the slider observe how the separation of Bragg disks and probe size change with convergence angle for graphene. 80 kV simulation.
```

For more complex analysis, there are two most common modes for experimental set-up: small and large convergence angle, which can be explored in [](#fig_stem_graphene). Using a small convergence angle leads to large separation between Bragg disks in reciprocal space, which is ideal for orientation and strain experiments. For orientation mapping, the diffraction pattern at each position is compared to a library of patterns from a reference structure to solve for rotation of the crystal. Strain mapping is accomplished by comparing to a reference diffraction pattern and computing  subtle shifts in Bragg disks resulting from structural disorder. The structure of amorphous materials can also be studied with this configuration by computing the radial distribution function. Careful choice of simulation parameters is imporant for this set-up. The small convergence angle creates a large real space probe size, which requires a large cell to hold the probe. At the same time, the maximum scattering angle needs to be large enough to capture enough Bragg disks for structural identification.

A large convergence angle creates a small probe in real space, which is ideal for atomic resolution imaging, and the overlap between diffracted disks transfers information for phase reconstruction experiments. There are a variety of ways to reconstruct the phase of the sample from 4D-STEM experiments, and the ideal approach will depend on the structure, field of view, and simulation parameters. The simpilest approach is a DPC reconstruction. As compared to the experimental set-up in the [DPC](#differential-phase-contrast) section, the segmented detectors are replaced with a pixelated detector for more accurate center of mass and thus phase determination.  More complicated reconstruction algorithms include parallax, where virtual images from each pixel in the central disk are aligned to create a high-signal to noise phase reconstruction. Interative electron ptychography referes to a family of algoirthms where the probe and object are reconstructed from the 4D-STEM dataset. In this method the maximum scattering angle is set by the highest spatial frequency recorded in the detector, so it is possible to form an image with higher resolution than real space sampling. More information about phase contrast reconstructions can be found at [].