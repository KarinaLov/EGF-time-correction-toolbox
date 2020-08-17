# EGF-time-correction-toolbox
Matlab toolbox for estimating the Green's function from cross correlating ambient seismic noise and measure time delays in seismic data.

The codes presented here were initially developed as a part of my master thesis "Measuring seismic station timing errors from ambient noise" (http://bora.uib.no/handle/1956/18942) where the main objective was to create a tool that estimates the Green's function from ambient seismic noise and use the Green's function to measure instrumental timing errors in seismic data. These codes have been rewritten to functions and are put together as toolbox. The master thesis and codes were written at the Department of Earth Science at the University of Bergen. The codes have also been used in https://pubs.geoscienceworld.org/ssa/srl/article/91/4/2247/586592/When-Clocks-Are-Not-Working-OBS-Time-Correction.

The toolbox contains two main functions, one for estimating the Green's function and one for measuring time shifts. In addition, sub- functions for preparing, processing and cross correlating the input data, and side functions for analysing, inverting and plotting the results, are included. The input values are all defined and can be changed in a settings text file, it is therefore only necessary to specify the name of the settings file as input to the main functions.

Most of the functions are written in MATLAB, but some are also written with SAC and Linux shell. Data samples are distributed together with the codes to demonstrate how the codes work. 

NB! The functions rdsac.m, mksac.m and ReadMSEEDFast.m reads mseed-files must be downloaded and placed in the toolbox before the functions can be used.  They can be found on these links: https://se.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac-read-and-write-sac-seismic-data-file and
https://se.mathworks.com/matlabcentral/fileexchange/46532-readmseedfast-filename

