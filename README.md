MTEX-texture-block-analysis
-----------

MTEX scripts for sequentially cropping electron backscatter diffraction (EBSD) maps, producing a grid of equivalent-sized squares (in x and y), for calculating and plotting the crystallographic texture variation across a sample.

Includes MTEX scripts to calculate the ODFs, calculate texture intensity values and plot pole figures for each individual square.

Contents
-----------
    
1. `texture_block_analysis.m` An MTEX script to segment the EBSD map into a grid of equivalent sized squares, to calculate crystallographic texture for the alpha phase, plot pole figures, plot ODFs and calculate texture strength values.

2. `texture_block_analysis_beta.m` An MTEX script to segment the EBSD map to calculate crystallographic texture for the beta phase.

MATLAB Setup
-----------

The scripts have been tested with MATLAB Version R2019b.

MTEX Installation
-----------

The MTEX toolbox can be downloaded from [here](https://mtex-toolbox.github.io/download), which also includes instructions for installing MTEX and troubleshooting any issues.

The scripts have been tested with MTEX Version 5.3, but should only require minor adjustments to work with future versions.

Some minor changes to the base MTEX code have been made to produce better looking figures. See [here](https://lightform-group.github.io/wiki/software_and_simulation/mtex-nice-figures) for more information.
