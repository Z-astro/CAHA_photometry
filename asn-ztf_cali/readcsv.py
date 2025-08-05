#! /usr/bin/env python8
# -*- coding: utf-8 -*-

import sys
import numpy as np
import os, glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


plt.rcParams["font.family"] = "serif"
mpl.rcParams["text.usetex"] = "True"
mpl.rcParams["font.size"] = 18
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["axes.labelpad"] = 10


def calmeanstd(mag, err):
    meanmag = sum(mag) / float(len(mag))
    mag_tmp = set(mag)
    if len(mag_tmp) > 1:
        sigma = float(np.std(mag, ddof=1))
        if sigma == 0:
            sigma = float(err[0])
    else:
        sigma = np.sqrt(np.sum(err**2) / len(err) ** 2)  # float(err[0])
    err_mean = np.sqrt(np.sum(err**2) / len(err) ** 2)
    if len(mag_tmp) > 1:
        newerr1 = np.sqrt(err_mean**2 + sigma**2 / len(err))
        newerr2 = np.sqrt(err_mean**2 + sigma**2)
    else:
        newerr1 = err_mean
        newerr2 = err_mean
    return meanmag, sigma, newerr1, newerr2


def combine_days_guo(jd, diffmag, differr, days):
    nbin = days
    data_dict = {"jd": jd, "mag": diffmag, "err": differr}
    combine_df = pd.DataFrame(data_dict)
    combine_df["cut_group"] = pd.cut(
        combine_df.jd, bins=np.arange(np.min(jd), np.max(jd) + nbin, nbin), right=False
    )
    grouped = combine_df.groupby(by="cut_group")
    new_jd = []
    new_diffmag = []
    new_sigma = []
    for tempname, tempgroup in grouped:
        if len(tempgroup) > 0:
            new_jd.append(np.sum(tempgroup["jd"].values) / float(len(tempgroup)))
            tmpdiffmag, tmpsigma, newerr1, newerr2 = calmeanstd(
                tempgroup["mag"].values, tempgroup["err"].values
            )
            new_diffmag.append(tmpdiffmag)
            new_sigma.append(tmpsigma / len(tempgroup) ** 0.5)
    return np.array(new_jd), np.array(new_diffmag), np.array(new_sigma)


def median_filter(x, n):
    """
    n = 3, 5, 7, 9...
    """
    nedge = (n - 1) // 2
    x1 = np.ones(nedge) * x[0]
    x2 = np.ones(nedge) * x[-1]
    newx = np.append(x1, x)
    newx = np.append(newx, x2)
    median = []
    for i in range(len(x)):
        y = np.median(newx[i : i + n])
        median.append(y)
    median = np.array(median)
    return median


def removeBadPoints(mjd, mag, err, crit=3, niter=20, n=9):
    for k in range(niter):
        jump = len(mag)
        mag_tmp = median_filter(mag, n)
        res = mag - mag_tmp
        flag = np.where(abs(res) <= crit * np.std(res))
        mjd, mag, err = mjd[flag], mag[flag], err[flag]
        if len(flag) == jump:
            break
    return mjd, mag, err


def ztf_lc(ztf_file, gap=1):
    df_ztf = pd.read_csv(ztf_file)
    df_ztf = df_ztf[df_ztf["catflags"] == 0]
    jd, mag, magerr = (
        df_ztf["mjd"].values,
        df_ztf["mag"].values,
        df_ztf["magerr"].values,
    )
    jd, mag, magerr = removeBadPoints(jd, mag, magerr)
    # jd_c, mag_c, magerr_c = combine_days_guo(jd, mag, magerr, gap)  # 按gap天数合并
    flux = 10 ** ((13.9 - mag) / 2.5)
    fluxerr = flux * magerr * np.log(10) / 2.5
    return [jd, flux, fluxerr]


def asas_lc(asas_file, Filter, gap=1):
    df_asas = pd.read_csv(asas_file)
    # print(df_asas["mag_err"].value_counts())  # 统计“mag_err”一列中有各值分别有多少个
    badPointsIndex = df_asas[
        (df_asas.mag_err == 99.990)
    ].index.tolist()  # 找出mag_err == 99.990的索引
    df_asas = df_asas.drop(
        index=badPointsIndex
    )  # 使用drop函数以及index参数删除上面找到的索引行
    df_asas = df_asas.reset_index(
        drop=True
    )  # 使用参数设drop=True删除旧序列，将标签重新从零开始顺序排序
    df_asas["mag"] = df_asas["mag"].astype(np.float64)
    df_asas["mag_err"] = df_asas["mag_err"].astype(
        np.float64
    )  # astype转换mag_err一列的格式
    # print(df_asas["Filter"].value_counts())  # 统计测光Filter(g/V)分别多少个点
    df_asas = df_asas[df_asas["Filter"] == Filter]  # 筛选出需要的测光Filter(g/V)
    df_asas = df_asas.reset_index(
        drop=True
    )  # 使用参数设drop=True删除旧序列，将标签重新从零开始顺序排序
    # print(df_asas)
    jd, mag, magerr = (
        df_asas["HJD"].values - 2400000.5,
        df_asas["mag"].values,
        df_asas["mag_err"].values,
    )
    jd, mag, magerr = removeBadPoints(jd, mag, magerr)
    # jd_c, mag_c, magerr_c = combine_days_guo(jd, mag, magerr, gap)  # 按gap天数合并
    # flux_c = 10 ** ((13.9 - mag_c) / 2.5)
    # fluxerr_c = flux_c * magerr_c * np.log(10) / 2.5
    flux = 10 ** ((13.9 - mag) / 2.5)
    fluxerr = flux * magerr * np.log(10) / 2.5
    return [jd, flux, fluxerr]


def write_txt_old(objnm):
    objpath = "/home/yao/astroWORK/CAHA_LC/AZ/" + objnm
    # [jd_zg, flux_zg, fluxerr_zg] = ztf_lc(objpath + 'zg.csv')
    [jd_zr, flux_zr, fluxerr_zr] = ztf_lc(objpath + "zr.csv")
    [jd_av, flux_av, fluxerr_av] = asas_lc(objpath + "asas.csv", "V")
    [jd_ag, flux_ag, fluxerr_ag] = asas_lc(objpath + "asas.csv", "g")
    # jd_zg, jd_zr, jd_av, jd_ag = jd_zg - 50000, jd_zr - 50000, jd_av - 50000, jd_ag - 50000
    jd_zr, jd_av, jd_ag = jd_zr - 50000, jd_av - 50000, jd_ag - 50000
    # jd_av, jd_ag = jd_av - 50000, jd_ag - 50000
    txtpath = "/home/yao/astroWORK/CAHA_LC/cali/" + objnm + ".txt"
    txt = open(txtpath, "wt")
    agformat = "# ASAS-SN-g %d\n" % jd_ag.shape[0]
    txt.write(agformat)
    for i in range(len(jd_ag)):
        rowformat = "%.6f   %.6f   %.6f\n" % (jd_ag[i], flux_ag[i], fluxerr_ag[i])
        txt.write(rowformat)
    zrformat = "# ZTF-r %d\n" % jd_zr.shape[0]
    txt.write(zrformat)
    for i in range(len(jd_zr)):
        rowformat = "%.6f   %.6f   %.6f\n" % (jd_zr[i], flux_zr[i], fluxerr_zr[i])
        txt.write(rowformat)
    # zgformat = '# ZTF-g %d\n' % jd_zg.shape[0]
    # txt.write(zgformat)
    # for i in range(len(jd_zg)):
    #     rowformat = '%.6f   %.6f   %.6f\n' % (jd_zg[i], flux_zg[i], fluxerr_zg[i])
    #     txt.write(rowformat)
    avformat = "# ASAS-SN-V %d\n" % jd_av.shape[0]
    txt.write(avformat)
    for i in range(len(jd_av)):
        rowformat = "%.6f   %.6f   %.6f\n" % (jd_av[i], flux_av[i], fluxerr_av[i])
        txt.write(rowformat)
    txt.close()

    # plot all lc
    fig = plt.figure(figsize=(16, 8))
    ax1 = fig.add_axes([0.1, 0.33, 0.8, 0.45])
    ax1.errorbar(
        jd_ag, flux_ag, yerr=fluxerr_ag, fmt=".", elinewidth=0.2, ms=6, c="blue"
    )
    ax1.errorbar(jd_av, flux_av, yerr=fluxerr_av, fmt=".", elinewidth=0.2, ms=6, c="k")
    # ax1.errorbar(jd_zg, flux_zg, yerr=fluxerr_zg, fmt='.', elinewidth=0.2, ms=6, c="purple")
    ax1.errorbar(
        jd_zr, flux_zr, yerr=fluxerr_zr, fmt=".", elinewidth=0.2, ms=6, c="red"
    )
    ax1.legend(loc=2, ncol=4, frameon=False, fontsize=8)
    ax1.minorticks_on()
    ax1.tick_params(right="on", left="on", which="both")
    ax1.set_ylabel(r"$\rm{Mag}$")
    ax1.set_xlabel("MJD-50000(days)")
    plt.show()


def write_txt(objnm):
    alltypes = ["asas", "zr", "zg", "zi"]
    asfilters = ["g", "V"]
    colordic = {"g": "k", "V": "blue", "zr": "red", "zg": "green", "zi": "orange"}
    # lcfileDir = '/home/yao/astroWORK/CAHA_LC/AZ/'
    lcfileDir = f"../../asn-ztf/obj_{objnm}"
    lcfilelst = sorted(glob.glob(os.sep.join((lcfileDir, objnm + "*.csv"))))
    # print(lcfilelst)
    # lctypelst = [
    #     lcfilelst[i][len(lcfileDir) + len(objnm) : -4] for i in range(len(lcfilelst))
    # ]
    lctypelst = [
        os.path.basename(filepath)[len(objnm):-4] for filepath in lcfilelst
    ]
    print(lctypelst)  # 这里默认ASAS-g是在第一位
    type1 = "asas"  # 指定以哪个type为基准，就把它放到第一位
    if type1 in lctypelst:
        index = lctypelst.index(type1)  # 找到type1元素的索引
        lctypelst = (
            [lctypelst[index]] + lctypelst[:index] + lctypelst[index + 1 :]
        )  # 将type1元素放到第一位，其他元素位置顺序不变
    print(lctypelst)
    # for i in range(len(alltypes)):
    #     if alltypes[i] not in lctypelst:
    #         print('<<<< No (%s) lc for (%s) >>>>' % (alltypes[i], objnm))

    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_axes([0.1, 0.25, 0.8, 0.6])
    # txtnm = "/home/yao/astroWORK/CAHA_LC/cali/" + objnm + ".txt"
    txtnm = f"../../asn-ztf/obj_{objnm}/{objnm}.txt"
    txt = open(txtnm, "wt")
    for i in range(len(lctypelst)):
        lcfile = os.sep.join((lcfileDir, objnm + lctypelst[i] + ".csv"))
        print(lcfile)
        if lctypelst[i] == "asas":
            for j in range(len(asfilters)):
                [jd, flux, fluxerr] = asas_lc(lcfile, asfilters[j])
                jd -= 50000
                ax.errorbar(
                    jd,
                    flux,
                    yerr=fluxerr,
                    fmt=".",
                    elinewidth=0.2,
                    ms=6,
                    c=colordic[asfilters[j]],
                    label="ASAS-SN-" + asfilters[j],
                )
                ax.legend(loc=2, ncol=4, frameon=False, fontsize=14)
                row0format = "# ASAS-SN-" + asfilters[j] + " %d\n" % jd.shape[0]
                txt.write(row0format)
                for k in range(jd.shape[0]):
                    rowformat = "%.6f   %.6f   %.6f\n" % (jd[k], flux[k], fluxerr[k])
                    txt.write(rowformat)
        else:
            [jd, flux, fluxerr] = ztf_lc(lcfile)
            jd -= 50000
            ax.errorbar(
                jd,
                flux,
                yerr=fluxerr,
                fmt=".",
                elinewidth=0.2,
                ms=6,
                c=colordic[lctypelst[i]],
                label="ZTF-" + lctypelst[i][-1],
            )
            ax.legend(loc=2, ncol=4, frameon=False, fontsize=14)
            row0format = "# ZTF-" + lctypelst[i][-1] + " %d\n" % jd.shape[0]
            txt.write(row0format)
            for k in range(jd.shape[0]):
                rowformat = "%.6f   %.6f   %.6f\n" % (jd[k], flux[k], fluxerr[k])
                txt.write(rowformat)
    txt.close()
    ax.minorticks_on()
    # fig.text(0.5, 0.045, r'$L_{5100}\ \rm(erg\ s^{-1})$', va='center', ha='center', fontsize=20)  # 设置整体的xlabel
    # fig.text(0.035, 0.5, r'$R_{\rm H\beta}\ \rm(lt-days)$', va='center', rotation='vertical', fontsize=20)  # 设置整体的ylabel
    ax.set_xlabel("MJD - 50000 (days)", fontsize=22)
    ax.set_ylabel(r"$\rm{Mag}$", fontsize=22)
    # plt.show() # open it only when you update single source. Or it will block the cali_targets.sh
    plt.close()


def main():
    objnm = sys.argv[1]
    # objnm = 'IC4329A'
    write_txt(objnm)


if __name__ == "__main__":
    main()
