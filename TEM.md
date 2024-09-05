---
title: Transmission Electron Microscopy
numbering:
  enumerator: 1.%s
---

(tem_sims)=
## TEM Simulations with Parallel Illumination

#todo: add missing citations

### Wavefunctions
After building an atomic potential as described in the [](#algorithms_page), the first step in a TEM simulation is to choose the wavefunction ({math}`\Psi`) for the simulation. The simplest case for the incident beam is to set {math}`\Psi` to unity everywhere in the plane, which means perfectly even illumination across the sample. However, it is possible to introduce more complications, such as a slightly [titled plane wave](https://abtem.readthedocs.io/en/main/user_guide/walkthrough/multislice.html#small-angle-beam-tilt). The sampling is set by the gridpoint and extent as described in seciton [](#sim_inputs_page).

### Imaging

A TEM image is simulated by propagating the wavefunction through the specimen potential using the multislice algorithm, which calculates how the wave evolves due to scattering by the specimen atoms and the propagation through it. The resulting exit-wave is complex, but can be visualized via its intensity. For a more realistic image, a [contrast transfer function](#CTF_page) can be applied to model the optics of the microscope; note that this is done after the computationally time-consuming multislice run which described the physics of hte interaction.

In [](#fig_tem_Au_potential_wave_image) we show an interactive visualization of the potential of gold with a lattice constant of 4.08 Å in the <100> zone axis, the corresponding exit wave function, and the resulting image as a function of depth through the specimen.

```{figure} #app:tem_Au_potential_wave_image
:name: fig_tem_Au_potential_wave_image
:placeholder: ./static/tem_Au_potential_wave_image.png
**Visualization of slicing through the specimen for the potential, the exit wave, and the image with a CTF applied.**
```

%#### Imaging example 
%
%```{figure} #app:tem_imaging
%:name: fig_tem_phase
%:placeholder: ./static/tem_imaging.png
%**TEM imaging of SrTiO$_3$ grains**: 
%```

### Electron diffraction patterns

Instead of an image, we can instead simulate a selected area diffraction (SAD) experiment by using the `DiffractionPatterns` method. We use `block_direct=True` to block the direct beam: it typically has a much higher intensity than the scattered beams, and thus it is typically not possible to show it on the same scale.

You may wonder; why do the diffraction spots look like squares? This is because the incoming wave function is a periodic and infinite plane wave, hence the intersection with the Ewald sphere is pointlike. However, since we are discretizing the wave function on a square grid (i.e. pixels), the spots can only be as small as single pixels. In real SAD experiments, the spot size is broadened due to the finite extent of the crystal as well instrumental effects.

We can use the `index_diffraction_spots` method to create a represention of SAD patterns as a mapping of Miller indices to the intensity of the corresponding reflections. The *conventional* unit cell have to be provided in order to index the pattern, we can provide this as the unit cell of the gold crystal we created earlier, we cannot use the the repeated cell.

In [](#fig_tem_Au_diffraction) we show an interactive visualization of the diffraction intensities and indexed diffraction spots as function of the depth through the specimen. We see that the {100} reflections are extinguished, as is expected from the selection rules of an F-centered crystal. We can also observe that the <220> spots end up with significantly higher intensity than the <200> spots; this is due to dynamical scattering — which is accounted for by the multislice algorithm.

```{figure} #app:tem_Au_potential_wave_image
:name: fig_tem_Au_diffraction
:placeholder: ./static/tem_Au_diffraction.png
**Visualization of redistribution of diffraction intensity as function of depth through an Au <100> specimen, and the Miller indexing of the resulting diffraction spots.**
```

%#### Diffraction example
%
%```{figure} #app:tem_diffraction
%:name: fig_tem_diffraction
%:placeholder: ./static/tem_diffraction.png
%**TEM diffraction of STO as a function of thickness**: 
%```

### Contrast transfer and phase contrast STEM
Thus far we have been considering how to form images and diffraction patterns with perfect incident illumination. However, often we're interested in seeing how aberrations or other beam modifications impacts imaging conditions. There are a variety of aberration functions ({math}`\chi(\bm{k})`) we may be interested in including as described in [](#CTF_page).

In TEM, aberrations modify the exit wave ({math}`\Psi_{exit}`) after the multislice simulation: 

```{math}
:label: eq:TEM_aberrations

\Psi_{image}(\bm{k}) = \Psi_{exit}(\bm{k}) \mathrm{e}^{-i\chi(\bm{k})},
```

One example of the importance of aberrations and wavfunction modification is shown in [](#tem_contrast). Although the wavfunction at the imaging plane ({math}`\Psi_{image}(\bm{k}`)) may be complex, only the intensity of the exit wave is measured by detectors. Weakly scattering samples, impart very little ampltidue contrast on the incident beam, so most of the structural information is encoded in the phase of the exit wave. This is especially well known in the case of imaging of biological materials, as these structures are beam-sensitive and composed of low atomic number (weakly scattering elements). This effect is evident in the left pannel of [](#tem_contrast), where as the dose decreases, it is very challenging to observe the simulated covid spike protein.

```{figure} #app:tem_contrast
:name: fig_tem_phase
:placeholder: ./static/tem_contrast.png
**TEM imaging of COVID spike protein**: (left) In focus, (middle) defocued, and (right) Zernike phase plate TEM images. [Covid structure](https://www.rcsb.org/structure/3jcl) from {cite}`walls2016cryo`
```

For the same, dose the image contrast can be increased by modfying the wavefunctions through the introduction of aberrations. This is most often done with defocus, which is illustrated in the middle pannel. The resulting image has improved contrast, making it easier to visualize the covid spike protein. The middle pannel shows the defocused image, but the aberrations can often be corrected for later in post-processing. 

The phase of the spike protein can also be image through the introduction of a phase plate after the sample. This simulation ([](#tem_contrast) right) shows a Zernike phase plate, where electrons scattered to higher angles, gain an extra phase shift. This high-contrast image is free of aberrations, making it easy to see the spike protein. However, changing the slides on the widget also highlgihts some of the challenges of phase contrast imaging with a phase plate. The best transfer of informaiton occurs with a phase plate that approximtes a delta function as closely as possible, meaning all the electrons exepct those in the very center of the aperature get an extra phase shift. However, as the size of the hole increases or the phase shift changes, the transfer of information decreaes. 
