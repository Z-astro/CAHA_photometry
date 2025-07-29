#!/usr/bin/env bash

# 重要！必须加载 conda 初始化
source ~/miniconda3/etc/profile.d/conda.sh

start_time=$(date +%s)  # 开始时间

curpath=`pwd` 
# 注意：以下路径是硬编码的，如果您的目录结构不同，请务必修改
rootdir=/mnt/WDsata/CAHA
baseDirForScriptSelf=/mnt/WDsata/CAHA/codes/step2

current_date=$(date +%Y%m%d)
echo "当前日期：$current_date"

# 注意：以下路径是硬编码的
lcpath=/mnt/WDsata/CAHA/lightcurves/lightcurve_${current_date}/
mkdir -p ${lcpath}

for dirname in `ls obj* -d`;do
  echo "<<<<<${dirname}>>>>>"
  cd ${dirname}
  
  # 为 demo1.py 使用 astroconda 环境 (Python 3.7.16)
  (
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate astroconda
    # python ${baseDirForScriptSelf}/demo1.py
    python ${baseDirForScriptSelf}/alipy_align.py
  )
  
  # 检查 alipy_out 目录是否存在，如果不存在则跳过后续步骤
  if [ ! -d "alipy_out" ]; then
      echo "目录 alipy_out 不存在，跳过 ${dirname}"
      cd ${curpath}
      continue
  fi
  
  cd alipy_out/
  
  # 检查是否有FITS文件，如果没有则跳过
  if ! ls *.fits > /dev/null 2>&1; then
      echo "alipy_out 目录中没有 FITS 文件，跳过 ${dirname}"
      cd ${curpath}
      continue
  fi
      
  ls *.fits > objlst
  
  # 为 ca_photutils_py3.py 使用 lc_CAHA 环境 (Python 3.13.0)
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
done

end_time=$(date +%s)  # 结束时间
runtime=$((end_time - start_time))  # 计算脚本运行耗时
echo "脚本运行时间为：$runtime 秒"