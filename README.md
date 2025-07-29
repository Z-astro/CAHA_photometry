# CAHA_photometry
The photomtery measurement pipeline for CAHA 2.2m Telescope (CAFOS) Reverberation Mapping Survey.

CAHA (Centro Astronómico Hispano en Andalucía) is located in Calar Alto Observatory in Spain. The SARM group in IHEP uses CAFOS 2.2m Telescope for RM survey. This pipeline is created for photometry of the obseravations.

The first man who governed the CAHA photometry was Zhi-Xiang Zhang, then was Wei-Jian Guo, then was Zhu-Heng Yao. Now it's Hao Zhang. (By Jul 2025)

The whole pipeline is dividied into 2 parts. Remember install the conda environment `astroconda` (python==3.7.16) which supports the `pyraf` in python3.

## Step 1
In this step, we read the original data from SFTP server in CAHA. We put the sources into their own folders and cut the images to the scale we want. Originally, the most important step is conducted by `reduce_img.py`. But this scripy is designed for python==2.7 which is too old. Fortunately it's not difficult to upgrade it for python3. As one can see, just assign the terms don't need a correction in function `trim` explicitly, the script can run. All the scripts in this step are expected to be executed in `astroconda`.

## Step 2
In this step, we run the aperture photomtery using `photutils` to obtain lightcurves for every source. First we have to align all the images to the reference image so that flux computation is available. Image alignment is completed using `alipy_align.py` in `astroconda`. The name was `demo1.py` but Hao Zhang changed it because it looks so unprofessional. After that, run `ca_photutils_py3.py` in an environment with python==3.13.0 or other versions that allows `photutils=2.2.0`. Last, run `FWHM_measurement.py` to calculate the FWHM of the target and Spec.Comparison star as the seeingo of the image.
