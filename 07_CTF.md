---
title: Wave Aberrations
numbering:
  enumerator: 1.%s
label : CTF_page
---

An ideal lens forms a spherical wave converging on or emerging from a single point. This means that a plane wave will be converted to a spherical wave converging at the focal point of the lens, and the image of a point source is also a point. In STEM we want the objective lens to produce the smallest possible probe and in HRTEM we want the objective lens to produce a perfect magnified image of the sample. In most cases, both of these objectives require that we minimize all aberrations as much as possible. 

However, while the last decades has seen enormous improvements in the optics of electron microscopes, they are far from ideal optical system. Imperfections causes the focused wave front to deviate from the ideal spherical surface. 

This deviation is typically expressed as a phase error or the aberration function, $\chi(\bm{k})$. Given a Fourier space wavefunction $\Psi_0(\bm{k})$ entering the lens, the wavefunction after passing through that lens can thus be expressed as 
```{math}
    \Psi(\bm{k}) = \Psi_0(\bm{k}) \mathrm{e} ^ {-i \chi(\bm{k})}.
```

Another, possibly more intuitive way of framing the phase error is the point spread function 
```{math}
    \mathrm{PSF}(\bm{r}) = F \mathrm{e} ^ {-i \chi(\bm{k})} ,
```
the point spread function describes how an imaging system responds to a point source. In STEM, the PSF would be the image of the probe, given an infinite objective aperture, in HRTEM the PSF would be how a perfect point source is imaged by an objective lens with an infinite collection angle.  

The phase error can be expressed in different ways, however, it is traditionally written as a series expansion. One popular choice is to write it as an expansion in polar coordinates
```{math}
    \chi(k, \phi) = \frac{2 \pi}{\lambda} \sum_{n,m} \frac{1}{n + 1} C_{n,m} (k \lambda)^{n+1} \cos\left[m (\phi - \phi_{n,m}) \right] ,
```
where $k$ is the Fourier space radial coordinate and $\phi$ is the corresponding azimuthal coordinate. The coefficient $C_{n,m}$ represents the magnitude of an aberration and $\phi_{n,m}$ gives a direction to that aberration. There are other ways of representing the aberration function, you might have seen the Cartesian representation where some parameters has an "a" and "b" version.

For an uncorrected microscope the dominant aberration is ther third order spherical aberration. Assuming the other non-symmetric components have been well aligned by the user, the contrast transfer function simplifies to:
```{math}
    \chi(k) \approx \frac{2\pi}{\lambda}\left( \frac{\lambda^2 k^2}{2} \Delta f + \frac{\lambda^4 k^4}{4} C_s \right) \quad .
```
Here we used the common aliases of the aberration coefficients, so $C10 = -\Delta f$ is the negative defocus and $C_{30} = C_s$ is the third order spherical aberration. 