This is a Matlab script for HiFi-SIM that allows for batch processing of raw SIM data inside a folder. It will process single-image tif files and stack files. It is compatible with raw SIM images as a series of planes (order: phases, angles) and as a 3x3 (phases, angles) mosaic (raw format from the Nikon N-SIM microscope). It can output the reconstructed widefiled images, classical Wiener reconstruction, and HiFi-SIM reconstruction.
To use, place the HiFiSIM_NC_batch.m and imgsave32seq.m files in the "Main_fun" folder of the HiFi-SIM Matlab code [available as a supplement](https://www.nature.com/articles/s41377-021-00513-w#Sec15) to the Hifi-SIM article:

Gang Wen, Simin Li, Linbo Wang, Xiaohu Chen, Zhenglong Sun, Yong Liang, Xin Jin, Yifan Xing, Yaming Jiu, Yuguo Tang & Hui Li.
High-fidelity structured illumination microscopy by point-spread-function engineering.
Light Sci Appl 10, 70 (2021).
[https://doi.org/10.1038/s41377-021-00513-w](https://www.nature.com/articles/s41377-021-00513)

HiFiSIM_NC_batch.m is the script file, and imgsave32seq.m is a modified saving routine that allows to control if images are appended to an existing stack or newly created. At the beginning of the HiFiSIM_NC_batch.m file you can:
- set the input folder path
- choose the outputs:
  - reconstruted widefield
  - Wiener reconstruction
  - HiFi-SIM reconstruction
- set various options for processing:
  - split slices into distinct output tifs,
  - estimate illumination pattern only for first file and use for other image files
  - estimate illumination pattern only for first slice and use for other slices within an image file
- set the optical configuration:
  - raw image pixel size in nm
  - objective NA
  - emission wavelength
- set the Hifi-SIM parameters:
    - attenuation stength (for optical sectioning)
    - damping factor
    - attenuation FWHM
- declare 3D-SIM for slice-by-slice processing (not tested, might be non-functional)
