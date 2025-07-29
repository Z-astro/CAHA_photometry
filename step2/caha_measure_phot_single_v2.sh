#!/usr/bin/env bash

# 重要！必须加载 conda 初始化
source ~/miniconda3/etc/profile.d/conda.sh

start_time=$(date +%s)  # 开始时间

curpath=`pwd` 
rootdir=/mnt/WDsata/CAHA
baseDirForScriptSelf=/mnt/WDsata/CAHA/codes/step2

current_date=$(date +%Y%m%d)
echo "当前日期：$current_date"

lcpath=/mnt/WDsata/CAHA/lightcurves/lightcurve_${current_date}/
mkdir -p ${lcpath}

# dirlst=(obj_IZw1)
# dirlst=(obj_Ark120 obj_SDSSJ082015 obj_SDSSJ094927)
# dirlst=(obj_SDSSJ082015)
# dirlst=(obj_SDSSJ094927)
# dirlst=(obj_CM01)
# dirlst=(obj_CM04)
# dirlst=(obj_IRAS13349)
# dirlst=(obj_Mrk509)
# dirlst=(obj_Mrk509B)
# dirlst=(obj_Mrk509R)
# dirlst=(obj_PDS456)
# dirlst=(obj_PG1126-041)
# dirlst=(obj_PG1700+518)
# dirlst=(obj_NGC4151)
# dirlist=(obj_PG1545+210)
# dirlist=(obj_PG2308+098)
# dirlist=(obj_PG0934+013)
# echo $dirlist
# dirname=${dirlist[*]}
dirname="obj_PG1435-067"

echo "<<<<<${dirname}>>>>>"
cd ${dirname}

# 在 astroconda 环境中运行 demo1.py. Now it's alipy_align.py
(
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate astroconda
    # python ${baseDirForScriptSelf}/demo1_py3.py
    # python ${baseDirForScriptSelf}/demo1.py
    python ${baseDirForScriptSelf}/alipy_align.py
)

# 检查 alipy_out 目录是否存在
if [ ! -d "alipy_out" ]; then
    echo "目录 alipy_out 不存在，跳过 ${dirname}"
    cd ${curpath}
    exit 1
fi

cd alipy_out/

# 检查是否有FITS文件
if ! ls *.fits > /dev/null 2>&1; then
    echo "alipy_out 目录中没有 FITS 文件，跳过 ${dirname}"
    cd ${curpath}
    exit 1
fi

ls *.fits > objlst

# 在 lc_CAHA 环境中运行 ca_photutils_py3.py
(
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate lc_CAHA
    python ${baseDirForScriptSelf}/ca_photutils_py3.py
    python ${baseDirForScriptSelf}/FWHM_measure.py
)

# 确保 lightcurve 目录存在再执行拷贝
if [ -d "../lightcurve" ]; then
    cp -r ../lightcurve/*obj.png ../lightcurve/*obj.lc ../lightcurve/*_fwhm_photutils.dat ../lightcurve/*_fwhm_photutils.png ${lcpath}
else
    echo "警告: 目录 ../lightcurve 未找到，无法拷贝结果文件。"
fi

cd ${curpath}

end_time=$(date +%s)  # 结束时间
runtime=$((end_time - start_time))  # 计算脚本运行耗时
echo "脚本运行时间为：$runtime 秒. ${dirname} 更新完成。"