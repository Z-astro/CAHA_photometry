import sys
import pycali


objnm = sys.argv[1]
# objnm = 'IZw1'
# txtnm = "/home/yao/astroWORK/CAHA_LC/cali/" + objnm + ".txt"
txtnm = f"../../asn-ztf/obj_{objnm}/{objnm}.txt"

cfg = pycali.Config()

cfg.setup(
    fcont=txtnm,  # fcont is a string
    fline=[],  # fline is a list, include multiple lines
    nmcmc=10000,
    ptol=0.1,
    scale_range_low=0.5,
    scale_range_up=2.0,
    shift_range_low=-1.0,
    shift_range_up=1.0,
    syserr_range_low=0.0,
    syserr_range_up=0.2,
    errscale_range_low=0.5,
    errscale_range_up=2.0,
    sigma_range_low=1.0e-4,
    sigma_range_up=1.0,
    tau_range_low=1.0,
    tau_range_up=1.0e4,
    fixed_scale=False,
    fixed_shift=False,
    fixed_syserr=True,
    fixed_error_scale=True,
    fixed_codes=[],  # fixed_codes is a list to specify the codes that need not to intercalibrate
    # e.g., [1, 3], will fix 1st and 3rd codes
)
cfg.print_cfg()

cali = pycali.Cali(cfg)  # create an instance
cali.mcmc()  # do mcmc
cali.get_best_params()  # calculate the best parameters
cali.output()  # print output
cali.recon()  # do reconstruction

# === plot results to PyCALI_results.pdf
# pycali.plot_results(cfg)

# === a simple plot
# pycali.simple_plot(cfg)
#
