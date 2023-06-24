---
title: Introduction
numbering:
  enumerator: 1.%s
---


% This is a direct citation from \cite{ophus2019four}, and this is an indirect citation \citep{ophus2019four}.

This article covers everything you need to know to simulate images for scanning/transmission electron microscopy (S/TEM) experiments, with a focus on a modern Python implementations. This text is aimed at both beginner and intermediate users, though experts may also learn some new concepts. The goal of the text is to enable readers to have an understanding of both the theory and practical applications of S/TEM image simulation. After reading, you should feel comfortable in devising, preparing and running S/TEM image simulations, and be able to critically assess simulated images in literature. To that end, this text contains interactive figures, with the hope that these may aide the reader in understanding/visualising the described concepts, as well as familiarising them with the associated code. While there are no formal pre-requisites and the text aims to be self-contained, a basic understanding of the Python programming language, linear algebra, and condensed matter physics will enhance the understanding of the text.

Indirect citation example:  {cite:p}`ophus2019four`

Direct citation example:  {cite:t}`ophus2019four`


The text can be broken into three sections: (1) Basic introduction to Python as well as crucial mathematical concepts. (2) Theory of S/TEM image simulation. (3) Practical implementations and considerations for image simulation. We finish by offer our perspective on the outlook for such simulations. The first section consists of a brief overview to the [Python programming language](#code_stuff), and instructions on how to run the interactive figures. This is followed by a primer on the [mathematical concepts](#math_stuff) that underpin image simulations. The second section addresses the theory of S/TEM image simulation, beginning with a discussion on the [underlying physics](#physics_stuff) of electron-atom interactions, after which the [algorithms](#electron-wavefunctions-ts) used to simulate images are detailed. In the final section the practicalities of image simulation are discussed, covering [creating simulation inputs](#sim_inputs), [TEM simulations](#tem_sims), [STEM simulations](#stem_sims), post-processing simulated images to create more experimentally realistic images, common errors and helpful tips. Finally we conclude with an outlook for image simulations.