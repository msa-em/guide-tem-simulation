---
title: Transmission Electron Microscopy
numbering:
  enumerator: 1.%s
---

(tem_sims)=
## TEM Simulations with Parallel Illumination

#todo: add missing citations

### Wavefunctions
After building an atomic potential as described in the [](#algorithms_page), the first step in a TEM simulation is to choose the wavefunction ({math}`\Psi`) for the simulation. The simpilest case for the incident beam is to set {math}`\Psi` to unity everywhere in the plane, which means perfectly even illumination across the sample. However, it is possible to introduce more complications, such as a slightly [titled plane wave](https://abtem.readthedocs.io/en/main/user_guide/walkthrough/multislice.html#small-angle-beam-tilt). The sampling is set by the gridpoint and extent as described in seciton [](#sim_inputs_page).

In [](#fig_potential_wave_image) we show an interactive visualization of the potential, the corresponding exit wave function, and the resulting image with a reasonable [contrast transfer function](#CTF_page) applied, as a function of the slices through the specimen.

```{figure} #app:tem_potential_wave_image
:name: fig_potential_wave_image
:placeholder: ./static/potential_wave_image.png
**Visualization of slicing through the specimen for the potential, the exit wave, and the image with a CTF applied.**
```

### Imaging 

```{figure} #app:tem_imaging
:name: fig_tem_phase
:placeholder: ./static/tem_imaging.png
**TEM imaging of SrTiO$_3$ grains**: 
```



### Diffraction 

```{figure} #app:tem_diffraction
:name: fig_tem_diffraction
:placeholder: ./static/tem_diffraction.png
**TEM diffraction of STO as a function of thickness**: 
```



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
