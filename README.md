# A practical guide for Simulation of Scanning and Transmission Electron Microscopy Experiments



### Authors

Colin Ophus, Alex Rakowski, Hamish Brown, Jacob Madsen, Toma Susi, et al.

### Abstract



### Information

This repository contains the source files for this manuscript. Manuscript will be submitted to Elements of Microscopy, published by the [Microscopy Society of America](https://www.microscopy.org). Manuscript platform powered by [https://curvenote.com/](Curvenote).




### License - CC-BY

<img src="https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by.png" alt="drawing" width="200"/>

CC-BY: This license allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, so long as attribution is given to the creator. The license allows for commercial use.

### Building Manuscript Using MyST/Curvenote

Follow the [installation guide](https://myst.tools/docs/mystjs/quickstart) and then type `myst start`. This will bring up a live preview of the LaTeX document.
Most LaTeX environments/macros are supported in the current LaTeX document and all notebooks.

- With the webserver running, change the document and it will reprocess the document (<200ms).
- If you come across any LaTeX macros that you would like supported, or are not working correctly. Raise an [issue here](https://github.com/executablebooks/mystjs/issues/new/choose).
- To add a new notebook, add it in the `notebooks` folder and add to the list in the [_toc.yml (table of contents)](./_toc.yml).

Note that the IPyWidgets do **not** currently render in the default Curvenote theme. These will be supported in the MSA theme.