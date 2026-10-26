---
title: Simulation Inputs
numbering:
  enumerator: 1.%s
label : sim_inputs_page
---

(specimen-models)=
## Specimen Models

This chapter introduces the Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) for creating specimen models for use in TEM image simulation.

ASE is a set of tools and Python modules for setting up, manipulating and visualizing atomic structures, which is used in conjunction with a large number of atomistic simulation codes, for example [GPAW](https://wiki.fysik.dtu.dk/gpaw/) for running DFT simulations. In this notebook, ASE is introduced in the context of running electron microscopy image simulations with [*ab*TEM](https://abtem.github.io/doc/intro.html).

###  The `Atoms` Object

The `Atoms` object defines a collection of atoms. To define `Atoms` from scratch, we need to specify at least three things:

* atomic positions,
* atomic numbers (or chemical symbols),
* a periodic cell.

For example, to create a basic model of the N<sub>2</sub> molecule, we could define:

```Python
atoms = ase.Atoms("N2", positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], cell=[6, 6, 6])
```

All these attributes of the `Atoms` object are stored in underlying NumPy arrays, which can be directly modified if desired. Convenient arithmetic operations such as `+` and `*` also directly work for the `Atoms` object, so structures can be easily extended and combined to create more complex specimens models.

#### Importing Structures from Files

ASE can import all common atomic-structure formats (full list [here](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)). Below we import a `.cif`-file defining a unit cell of strontium titanate (SrTiO<sub>3</sub>) that we provide with this text and will use in further examples.

```Python
srtio3 = ase.io.read("srtio3.cif")
```

### Manipulating Atoms
*ab*TEM always assumes that the imaging electrons propagate along the $z$-axis in the direction from _negative to positive_ coordinate values. Hence, to choose the zone axis, we need to manipulate the atoms so they are properly aligned.

ASE has many tools for manipulating structures, but one particularly useful one is the `surface` function, which can be used for creating a periodic surface (aligned with the $z$-axis) for a given set of Miller indices.

In the widget below, we have oriented the strontium titanate structure along the (110)-direction and created supercells out of it, with 2 Å of vacuum added at the top and bottom surfaces.

```{figure} #app:sto_supercell
:name: fig_sto_supercell
:placeholder: ./static/sto_supercell.png
**Interactive widget showing supercell construction for the STO(110) supercell**:
```

Since the positions and atomic numbers are just `NumPy` arrays, they can be modified in-place. Below, we create an SrTiO<sub>3</sub>/LaTiO<sub>3</sub> interface by changing the atomic numbers of the Sr atoms with a $y$-coordinate less than $7.5 \ \mathrm{Å}$ in a (3,4,10) supercell oriented along the (110) zone axis. This interface created from a will be later used for [STEM image simulations](#stem-image-simulation).

```python
sto_lto = repeated_srtio3.copy()
mask = sto_lto.symbols == "Sr"
mask = mask * (sto_lto.positions[:, 1] < 7.5)
sto_lto.numbers[mask] = 57
```

### Non-Orthogonal or Tilted Specimen
Most multislice implementations also require the supercell axes to be orthogonal, which can sometimes be non-trivial to achieve for a given crystal structure, but especially tricky for rotated crystals. In *ab*TEM, there are some tools to help achieve this as described in its  [documentation](https://abtem.github.io/doc/user_guide/tutorials/advanced_atomic_models.html#orthogonal-and-periodic-supercells). While [small beam tilts](#tem-wavefunctions) can be accommodated within the multislice formalism, for arbitrary rotations of a crystal or large-angle beam tilts, the [Bloch Wave method](#blochwave-method) may instead be more appropriate (as implemented in [py4DSTEM](https://github.com/py4dstem/py4DSTEM_tutorials/blob/main/notebooks/diffraction_02_dynamical.ipynb) and recently also [*ab*TEM](https://abtem.github.io/doc/user_guide/tutorials/blochwave.html)).

## Sampling

In any numerical implementation, continuous physical quantities such as potentials or wavefunctions have to be described on numerical grids. In *ab*TEM, these are represented on a rectangular grid of $N_x \times N_y$ grid points (`gpts`) or pixels. 

Given an orthogonal cell with the sidelengths $L_x$ and $L_y$, in the $x$ and $y$-direction, the real-space sampling (or dimensions of the pixels) is $\Delta x = L_x / N_x$ in $x$ and $\Delta y=L_y / N_y$.
The real-space coordinates thus take on only discrete values of

$$
    x_i = i \Delta x \quad i = 0,1, \ldots , N_x - 1, \\
    y_j = j \Delta y \quad j = 0,1, \ldots , N_y - 1. \\ 
$$

The Fourier transform of this grid of values (for example, a simulated exit wave) will also have $N_x \times N_y$ grid points. However, in reciprocal space (for example, a simulated diffraction pattern), the sampling is instead determined by the supercell dimensions (real-space extent of the potential) given by the inverse relations $\Delta k_x = 1/L_x$ and $\Delta k_y = 1/L_y$. The reciprocal space coordinates thus take on values

$$
    k_{x,i} = - k_{x,\mathrm{max}} + i \Delta k_x, \\
    k_{x,j} = - k_{y,\mathrm{max}} + j \Delta k_y,
$$

where the maximum spatial frequencies (reciprocal-space extent) are imposed by the real-space sampling

$$ 
    k_{x,\mathrm{max}} = \frac{1}{2\Delta_x} \quad \mathrm{and} \quad k_{y,\mathrm{max}} = \frac{1}{2\Delta_y} \quad .
$$

Finally, we want to stress the perhaps non-intuitive fact that although some codes allow the sampling of diffraction patterns to be separately set, due to the reciprocal pixel-wise numerical correspondence between the real and reciprocal spaces, the only true way to improve reciprocal-space sampling, i.e. decrease $\Delta k$, is to increase the size of the supercell in $x$ and $y$.

In *ab*TEM we choose to retain the direct correspondence between the real-space extent of the potential and the reciprocal-space sampling, and only use interpolation to *decrease* the number of pixels along one direction to reach [uniform sampling](https://abtem.github.io/doc/reference/api/_autosummary/abtem.measurements.DiffractionPatterns.html#abtem.measurements.DiffractionPatterns.interpolate). You can find some more information in the [*ab*TEM documentation](https://abtem.github.io/doc/user_guide/appendix/antialiasing.html).

(choosing-sampling)=
### Choosing the Grid Size and Sampling

The real-space sampling and the size of the simulation cell are the two most important numerical parameters in a simulation, and choosing them poorly can produce errors which are easy to miss {cite:p}`kirkland2020`.

The real-space sampling $\Delta x$ sets the highest scattering angle included in the simulation. 
Because the transmission operator multiplies the wavefunction by the transmission function in real space, which is a convolution in Fourier space, the multislice method can produce aliasing artifacts at the edges of Fourier space. 
To prevent these, most codes including *ab*TEM apply an antialiasing aperture which removes all spatial frequencies above $2/3$ of $k_{\mathrm{max}}$, giving a maximum scattering angle of
```{math}
:label: eq:max_angle
\alpha_{\mathrm{max}} = \frac{2}{3} \frac{\lambda}{2 \Delta x} = \frac{\lambda}{3 \Delta x}.
```
For example, at 200 kV ($\lambda \approx 0.0251$ Å) a sampling of $\Delta x = 0.05$ Å gives $\alpha_{\mathrm{max}} \approx 168$ mrad. 
Any electrons scattered beyond this angle are removed from the simulation. 
For STEM simulations, $\alpha_{\mathrm{max}}$ must therefore be larger than the outer angle of every detector we want to simulate, otherwise the simulated ADF signals will be too low. 
The sampling must also be fine enough to describe the sharply peaked atomic potentials, which is why typical values are between 0.02 and 0.05 Å.

The size of the simulation cell $L$ sets the reciprocal-space sampling $\Delta k = 1/L$. 
In STEM, the cell must be large enough to contain the full probe, including the probe tails and any spreading of the probe through the thickness of the specimen and from defocus. 
Because the multislice method uses periodic boundary conditions, a probe which is too large for the cell will overlap with its own periodic images, producing wrap-around artifacts. 
The cell must also be large enough to sample the diffraction pattern finely, since the number of pixels across the diameter of the bright field disk is $2 \alpha_{\mathrm{cutoff}} L / \lambda$, where $\alpha_{\mathrm{cutoff}}$ is the convergence semiangle. 
For example, a $20$ mrad probe at 200 kV in a $20$ Å cell gives approximately 32 pixels across the bright field disk, which is sufficient for most detector geometries. 
Small or finely structured detectors, such as those used for DPC or 4D-STEM, require finer sampling of the diffraction pattern and therefore larger cells.

The safest way to choose these parameters is to perform a convergence test: repeat a small test simulation with finer sampling and a larger cell, and confirm that the quantity of interest (for example an image contrast or an integrated detector intensity) no longer changes by more than the required precision. 
Note that the real-space sampling of the wavefunction is independent of the STEM probe step size, which we discuss in [](#imaging).

(frozen-phonons)=
## Thermal Vibrations and Frozen Phonons

Atoms in a real specimen vibrate about their equilibrium positions, even at cryogenic temperatures due to zero-point motion. 
A fast electron passes through a typical TEM specimen in approximately $10^{-16}$ s, which is much shorter than the period of atomic vibrations (approximately $10^{-13}$ s). 
Each electron therefore sees the atoms frozen in a single set of displaced positions. 
In the frozen phonon method, we perform many multislice simulations, each with a different random set of atomic displacements, and average the resulting intensities incoherently {cite:p}`loane_thermal_1991`.

In the simplest model, the Einstein model, each atom is displaced independently by a random Gaussian displacement along each direction, with a root-mean-square displacement $u$ for each element. 
These values can be taken from the Debye-Waller factors $B = 8 \pi^2 u^2$ measured by X-ray or neutron diffraction, and are typically on the order of $0.1$ Å at room temperature. 
More accurate models use correlated displacements from molecular dynamics or phonon calculations {cite:p}`forbes_thermal_2011`.

Averaging the intensities over many configurations produces thermal diffuse scattering (TDS), the diffuse intensity between and beyond the Bragg disks, which dominates the signal at high scattering angles. 
TDS is therefore essential for quantitative simulations of HAADF-STEM images, and also reduces the intensity of the Bragg beams and changes the contrast of TEM images and CBED patterns. 
Even though phonon scattering is an inelastic process, the energy lost by the fast electron is only tens to hundreds of meV, and the frozen phonon method describes the resulting intensity accurately {cite:p}`van_dyck_frozen_2009`.

Because the frozen phonon result is a Monte Carlo average, its statistical error decreases with the number of configurations $N$ as $1/\sqrt{N}$. 
The number of configurations needed depends on the measurement. 
Signals which integrate over large regions of the diffraction pattern, such as ADF images, typically converge with fewer configurations than individual diffraction patterns, where the TDS between the Bragg disks remains speckled until many configurations have been averaged. 
Larger simulation cells also require fewer configurations, because each configuration already contains more independently displaced atoms. 
As for the sampling, a practical approach is to increase $N$ (for example from 4 to 8 to 16 configurations) until the quantity of interest changes by less than the required precision. 
In *ab*TEM, frozen phonons are defined by wrapping the atomic model, for example

```python
frozen_phonons = abtem.FrozenPhonons(atoms, num_configs=16, sigmas=0.1)
```

where `sigmas` is the root-mean-square displacement $u$ in Å, which can also be given separately for each element.
