These scripts are demos for FWHM measurements using IRAF. The core scripy is `batch_psfmeasure_ori.py`. The pipeline is:

1. batch_psfmeasurement_ori.py
2. txt2csv_psf.py
3. merge_psf_csv.py
4. fwhm_plot.py
5. cut_target_speccomp.py

The test source is PHL1092. If you want others please change the source name. The wierd thing is, the `batch_psfmeasure_ori.py` runs only if you press 'Q'. So you must press 'Q' continuously if you want it proceeds. I believe the wisdom of descendants. I give up. I use `photutils.psf.fit_fwhm()`. 
