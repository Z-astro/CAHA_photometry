# CAHA_photometry
The photomtery measurement pipeline for CAHA 2.2m Telescope (CAFOS) Reverberation Mapping Survey.

CAHA (Centro Astronómico Hispano en Andalucía) is located in Calar Alto Observatory in Spain. The SARM group in IHEP uses CAFOS 2.2m Telescope for RM survey. This pipeline is created for photometry of the obseravations. The framework of this version is produced by Zhu-Heng Yao (2019--2024, PhD).

The first man who governed the CAHA photometry was Zhi-Xiang Zhang, then was Wei-Jian Guo, then was Zhu-Heng Yao. Now it's Hao Zhang. (By Jul 2025)

The whole pipeline is dividied into 2 parts. Remember install the conda environment `astroconda` (python==3.7.16) which supports the `pyraf` in python3. The page is here: https://astroconda.readthedocs.io/en/latest/getting_started.html . You may follow the install instructions on Zhihu: https://zhuanlan.zhihu.com/p/545324660 .

## Step 1
In this step, we read the original data from SFTP server in CAHA. We put the sources into their own folders and cut the images to the scale we want. Originally, the most important step is conducted by `reduce_img.py`. But this scripy is designed for python==2.7 which is too old. Fortunately it's not difficult to upgrade it for python3. As one can see, just assign the terms don't need a correction in function `trim` explicitly, the script can run. All the scripts in this step are expected to be executed in `astroconda`.

## Step 2
In this step, we run the aperture photomtery using `photutils` to obtain lightcurves for every source. First we have to align all the images to the reference image so that flux computation is available. Image alignment is completed using `alipy_align.py` in `astroconda`. The name was `demo1.py` but Hao Zhang changed it because it looks so unprofessional. After that, run `ca_photutils_py3.py` in an environment with python==3.13.0 or other versions that allows `photutils=2.2.0`. Last, run `FWHM_measurement.py` in the same environment as `ca_photutils_py3.py` to calculate the FWHM of the target and Spec.Comparison star as the seeingo of the image. This script uses `photutils.psf.fit_fwhm()` in version 2.2.0. The other way is using `pyraf`, but I can't make it run automically.

`alipy_align.py` uses `astroconda` because it still needs `pyraf`. The alignment depends on `alipy` package: https://github.com/japs/alipy . This package is compatible with python3. But it needs other 2 packages installed in `astroconda`:
1. `asciidata`: https://github.com/japs/astroasciidata .
2. `SExtractor`: https://sextractor.readthedocs.io/en/latest/Introduction.html . This is a command-line tool.
3. `f2n`: https://github.com/zhang-zhixiang/f2n . This package is used for transformation between .fits and .png. This is an original version for python==2.7. So Hao Zhang corrects it to make it compatible in `astroconda`. 'Gemini is God', Hao Zhang said. The package corrected is in the repository.
4. `pyfits`: https://pypi.org/project/pyfits/3.5/ . Although this package is merged with `astropy`, `alipy` needs it. Install version 3.5.

`ca_photutils_py3.py` and `FWHM_measurement.py` demands a modern python environment to make sure `photutils==2.2.0` is available. By Jul 2025, Hao Zhang uses python==3.13.0. What's more, `FWHM_measurement.py` has another version `FWHM_measurement_bkg2D.py`, this version utilizes different background subtraction method. `FWHM_measurement.py`'s background subtraction keeps the same as `ca_photutils_py3.py`. In most cases, we use `FWHM_measurement.py` because using the same method with aperture photometry seems rational. If the evil referees ask you use another version, then you go.

## File Tree Structure
All files are located in the folder:`CAHA/`. 

The codes in this repository is in folder:`CAHA/codes/`. 

Original data is in `CAHA/data/tmp/`. 

The individual folder for every source is in `CAHA/objects/obj_{name}`. In every folder, original images after cut are in `CAHA/objects/obj_{name}/images/`. Images after alignment are in `CAHA/objects/obj_{name}/alipy_out/`. Lightcurves and FWHM curves are in  `CAHA/objects/obj_{name}/lightcurve/`. `CAHA/objects/obj_{name}/aper.txt` is the information for background subtraction in aperture photometry. `CAHA/objects/obj_{name}/{name}.pos` is a text file records the comparision stars'position. Note that the target is the last one, the Spec.comp is the last two. Others are photometry comparison. `CAHA/objects/obj_{name}/updatelst` records the progress.

The result folder is `CAHA/lightcurves`. The daily handling results are here. The name structure is: `lightcurve_{date}`.

`objradec.lst` is the file contains all the sources' position information. If you want to add new sources, remember to update it.

## ASAS-SN/ZTF lightcurve calibration
Corresponding folder: `ans-ztf_cali/`. The pipeline is `cali_targets.sh`. Usually we don't have to update once for all. Mr.Hu will tell you some targets are required to be updated. In Yao's version, you must divide the objects into several parts, a group with all ZTF-bands available, a group with only zg&zr bands available...... In this version, the script enters every source's directory to execute the mission, so you just need to give a total objlist in the shell script. If there is no demand, this section will be updated every 6 months.

ASAS-SN lightcurve download: https://asas-sn.osu.edu/

ZTF lightcurve download: https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mission=irsa&submit=Select&projshort=ZTF

**Caution! DON'T USE SKY-PATROL API KEY TO DOWNLOAD ASAS-SN DATA! Just download it from the website.**

Tools: PyCALI. https://pycali.readthedocs.io/en/latest/ . Hao Zhang said he still didn't know why this package is so magic. Ask Yan-Rong Li when you meet bugs during the installation.

Please correct the file name to the form like: '{objname}_asas.csv' and '{objname}_zg/zr/zi.csv'. Data file is located in: /mnt/WDsata/asn-ztf/obj _{name}/



