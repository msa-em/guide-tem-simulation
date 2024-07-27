---
title: Post Processing
numbering:
  enumerator: 1.%s
---
STEM simulations usually requires some post-processing, we apply some of the most common steps post-processing step in this tutorial.

Interpolation
We can save a lot of computational effort by scanning at the Nyquist frequency [insert reference], but the result is quite pixelated. To address this, we can interpolate the images to a sampling of Å. abTEM’s default interpolation algorithm is Fourier-space padding, but spline interpolation is also available, which is more appropriate if the image in non-periodic.



