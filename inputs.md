---
title: Simulation Inputs
numbering:
  enumerator: 1.%s
label : sim_inputs_page
---

text

(specimen-models)=
## Specimen models

This chapter introduces the Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) for creating specimen models for use in TEM image simulation.

ASE is a set of tools and Python modules for setting up, manipulating and visualizing atomic structures, which is used in conjunction with a large number of atomistic simulation codes, for example [GPAW](https://wiki.fysik.dtu.dk/gpaw/) for running DFT simulations. In this notebook, ASE is introduced in the context of running electron microscopy image simulations with [*ab*TEM](https://abtem.github.io/doc/intro.html).

###  The `Atoms` object

The `Atoms` object defines a collection of atoms. To define `Atoms` from scratch, we need to specify at least three things:

* atomic positions,
* atomic numbers (or chemical symbols),
* a periodic cell.

For example, to create a basic model of the N<sub>2</sub> molecule, we could define:

`atoms = ase.Atoms("N2", positions=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], cell=[6, 6, 6])`

All these attributes of the `Atoms` object are stored in underlying NumPy arrays,   which can be directly modified if desired. Convenient arithmetic operations also directly work for the `Atoms` object, so structures can be easily combined to create more complex specimens.

#### Importing structures from files

ASE can import all common atomic-structure formats (full list [here](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)). Below we import a `.cif`-file defining a unit cell of strontium titanate (SrTiO<sub>3</sub>) that we provide with this text and will use in further examples.

`srtio3 = ase.io.read("srtio3.cif")`

### Manipulating atoms
*ab*TEM always assumes that the imaging electrons propagate along the $z$-axis in the direction from _negative to positive_ coordinate values. Hence, to choose the zone axis, we need to manipulate the atoms so they are properly aligned.

ASE has many tools for manipulating structures, but one particularly useful one is the `surface` function, which can be used for creating a periodic surface (aligned with the $z$-axis) for a given set of Miller indices.

In the widget below, we have oriented the strontium titanate structure along the (110)-direction, and interactively create supercells out of it, with 2 Å of vacuum added at the top and bottom surfaces.

```{figure} #app:sto_supercell
:name: fig_sto_supercell
:placeholder: ./static/sto_supercell.png
**Interactive widget showing supercell construction for the STO(110) supercell.
```