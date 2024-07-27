---
title: Post-Processing
numbering:
  enumerator: 1.%s
---

(stem_post)=
## STEM Post-Processing
STEM simulations usually requires some post-processing, we apply some of the most common steps post-processing step in this tutorial.

Scanning imaging modes such as STEM works by rastering an electron probe across a sample pixel by pixel and recording the scattering signal. The computational cost of the simulation is directly proportional to the number of scan pixels, each requiring a separate multislice simulation.

For periodic speciments, even though the potential needs to be large enough to fit the probe, there is no need to scan over repated unit cells as tiliing afterwards can yield the same result. For the post-processing examples below, we use an STO-LTO heterointerface as a specimen, as shown below in [](#fig_stem_specimen).

```{figure} #app:stem_specimen
:name: fig_stem_specimen
:placeholder: ./static/stem_specimen.png
A SrTiO<sub>3</sub>/LaTiO<sub>3</sub> (STO/LTO) interface model.
```

### Interpolation
We can save a lot of computational effort by scanning at the Nyquist frequency [insert reference], but the result is quite pixelated. To address this, we can interpolate the images to a sampling of 0.05 Å. *ab*TEM’s default interpolation algorithm is Fourier-space padding, but spline interpolation is also available, which is more appropriate if the image in non-periodic.

### Blurring
A finite Gaussian-shaped source will result in a blurring of the image. Vibrations and other instabilities may further contribute to the blur. We apply a Gaussian blur with a standard deviation of $0.35 \ \mathrm{Å}$ (corresponding to a source of approximately that size). Note that correctly including spatial and temporal incoherence is a bit more complicated and may be necessary for quantitative comparisons with experiment.

### Noise
Simulations correspond to the limit of infinite electron dose. We can emulate finite dose by drawing random numbers from a Poisson distribution for every pixel. We apply Poisson noise corresponding a dose per area of $10^5 \ \mathrm{e}^- / \mathrm{Å}^2$.

The different STEM post-processing steps can be explored in [](#fig_stem_processing).

```{figure} #app:stem_processing
:name: fig_stem_processing
:placeholder: ./static/stem_processing.png
Using the slider observe how different post-processing steps affect the scanned bright-field, medium-angle, and high-angle annular dark-field images of an STO/LTO heterostructure.
```