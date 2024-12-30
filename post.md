---
title: Post-Processing
numbering:
  enumerator: 1.%s
---

(stem_post)=
## STEM Post-Processing
STEM simulations usually require some post-processing. We apply some of the most common steps post-processing step in this tutorial.

For these examples, we use an STO/LTO heterointerface as a specimen. The structure was built earlier in the [simulation inputs](./sim_inputs.md) chapter, and simple BF/ADF images simulated in the chapter on [STEM](./STEM.md).

#### Interpolation
We can save a great deal of computational effort by scanning at the [Nyquist_frequency](wiki:Nyquist_frequency), which is information-theoretically guaranteed to be sufficient — but the result is visually quite pixelated. To address this, we can interpolate the images to a sampling of 0.05 $\mathrm{\AA}$. *ab*TEM’s default interpolation algorithm is Fourier-space padding, but spline interpolation is also available, which is more appropriate if the image in non-periodic.

#### Blurring
Standard multislice simulations are too idealized to describe a realistic experimental image. For example, a finite Gaussian-shaped source will result in a blurring of the image, and vibrations and other instabilities may further contribute to the blur. It is typical and convenient to approximate these by applying a Gaussian blur with a standard deviation of $0.35 \ \mathrm{\AA}$ (corresponding to a source of approximately that size). However, note that correctly including spatial and temporal incoherence is a bit more complicated and may be necessary for quantitative comparisons with experiment.

#### Noise
Analogous to the discussion in [](#id-tem-phase), STEM simulations are initially performed at infinite dose, and we need to add Poisson nose to reach more realistic conditions. In this case we add a dose per area of $10^5 \ \mathrm{e}^- / \mathrm{\AA}^2$ to form a more realistic image.

The different STEM post-processing steps can be explored in [](#fig_stem_processing).

```{figure} #app:stem_processing
:name: fig_stem_processing
:placeholder: ./static/stem_processing.png
Using the slider observe how different post-processing steps affect the scanned bright-field, medium-angle, and high-angle annular dark-field images of an STO/LTO heterostructure.
```