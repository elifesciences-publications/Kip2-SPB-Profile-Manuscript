Info
=====
Raw data and scripts used for data analysis and modelling in [[1]].

License
-------
Code and Data are provided under the [MIT License](LICENSE).

Authors
-------

    Lukas Widmer <lukas.widmer@bsse.ethz.ch> (Scripts, Data Analysis)
    Xiuzhen Chen <xiuzhen.chen@bc.biol.ethz.ch> (Image Data)


Installation & Usage
====================

Windows
-------
1. Install MATLAB (tested on R2016b-R2018a)
2. Install Git via TortoiseGit (https://tortoisegit.org/)
3. Git clone this repo

Mac OS X
--------
1. Install MATLAB
2. Install Git through package manager: e.g., `brew install git`
3. Git clone this repo

Linux (instructions here are for Ubuntu / Debian)
-------------------------------------------------
1. Install MATLAB
2. Install Git through package manager: `sudo apt-get install git`
3. Git clone this repo

Usage
=====
* Run [runDataAnalysis.m](src/runDataAnalysis.m) in MATLAB to perform the analysis of the experimental data, and generate the data and model figures.
* Run [runSampling.m](src/modelFitting/runSampling.m) with one of the modes specified therein to simulate the stochastic model and compute the measurement model. Warning: this is computationally intense, using a compute cluster is highly recommended.

References
==========
[1]: http://google.com
1. Xiuzhen Chen†, Lukas A. Widmer†, Marcel M. Stangier, Michel O. Steinmetz, Jörg Stelling*, Yves Barral*. (2019)  
Remote control of microtubule plus-end dynamics and function from the minus-end
† contributed equally  
\* corresponding authors: <joerg.stelling@bsse.ethz.ch>, <yves.barral@bc.biol.ethz.ch>  
[Submitted.]