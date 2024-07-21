---
title: Introduction
numbering:
  enumerator: 1.%s
---

This article covers everything you need to know to simulate images for [transmission electron microscopy](wiki:Transmission_electron_microscopy) (TEM) and [scanning TEM](wiki:Scanning_transmission_electron_microscopy) (STEM) experiments [@pennycook2011scanning; @carter2016transmission]. We focus on modern [Python](wiki:Python_(programming_language)) implementations centered around the [numpy](wiki:NumPy) {cite:p}`harris2020array` and [cupy](wiki:CuPy) {cite:p}`nishino2017cupy` libraries, for CPU and GPU calculations respectively. 


This text is aimed at both beginner and intermediate users, though experts may also learn some new concepts. The goal of the text is to enable readers to have an understanding of both the theory and practical applications of S/TEM image simulation. After reading, you should feel comfortable in devising, preparing and running S/TEM image simulations, and be able to critically assess simulated images in literature. To that end, this text contains interactive figures, with the hope that these may aide the reader in understanding/visualising the described concepts, as well as familiarising them with the associated code. While there are no formal pre-requisites and the text aims to be self-contained, a basic understanding of the Python programming language, linear algebra, and condensed matter physics will enhance the understanding of the text.


The text can be broken into three sections:
1. Basic introduction to Python as well as crucial mathematical concepts.
2. Theory of S/TEM simulations.
3. Practical implementations and considerations for simulation. 

We finish by offer our perspective on the outlook for such simulations.
The first section consists of a brief overview to the [Python programming language](./02_code.md), and instructions on how to run the interactive figures.
This is followed by a primer on the [mathematical concepts](./03_math.md) that underpin image simulations. 
The second section addresses the theory of S/TEM image simulation, beginning with a discussion on the [underlying physics](./04_physics.md) of electron-atom interactions, after which the [algorithms](./05_algorithms.ms) used to simulate images are detailed. 
In the final section the practicalities of image simulation are discussed, covering [creating simulation inputs](./06_sim_inputs.md), [TEM simulations](./07_TEM.md), [STEM simulations](./08_STEM.md), post-processing simulated images to create more experimentally realistic images, common errors and helpful tips.
Finally we conclude with an outlook for image simulations.
