#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# run in the directory: /mnt/WDsata/CAHA/codes/asn-ztf_cali
start_time=$(date +%s)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate astroconda
echo 'activate conda env: astroconda (python==3.7.16)'

# pwdDir='pwd'
# baseDirForScriptSelf=$(cd "$(dirname "$0")"; pwd)  # 当前shell脚本所在目录
baseDirForScriptSelf=/mnt/WDsata/CAHA/codes/asn-ztf_cali # 硬编码路径，注意迁移时需要修改

# 2025.8.2 CAHA planlist

# ZTF with zg zr zi:
# IRAS13349 Mrk6 Mrk335 Mrk509 Mrk1239 NGC4151 OI371 PG0003+158 PG0043+039 PG0804+761 PG0844+349 PG0934+013 PG1011-040 PG1103-006 PG1126-041 PG1226+023 PG1244+026 PG1259+593
# PG1307+085 PG1354+213 PG1411+442 PG1416-129 PG1543+489 PG1700+518 PG2112+059 PG2130+099 PHL1092 PHL1811 SDSSJ094927 SDSSJ224028

# ZTF with zg zr (no zi):
# Ark120 IZw1 NED02 PDS456 PG1435-067 PG1545+210 PG1612+261 PG1704+608 SDSSJ082015 

# ZTF with zi zr (no zg):
# NGC7603 

# no ZTF data:
# PG2308+098
objlst=(PG1226+023 PHL1092 PG1416-129 PG1354+213 PG1307+085 PG1103-006 PG0934+013 PG1011-040 PG1411+442 PG1259+593 PG1244+026 
PG1612+261 PG1545+210 PG1704+608 SDSSJ082015)

current_date=$(date +%Y%m%d)
echo "当前日期：$current_date"

lcpath=/mnt/WDsata/CAHA/asn-ztf/lc_cali_update/lc_cali_${current_date}/ # 硬编码路径，注意迁移时需要修改
mkdir -p ${lcpath}
if [ $? -eq 0 ]; then
    echo "Successfully created directory: ${lcpath}"
else
    echo "Error creating directory: ${lcpath}" >&2
    exit 1
fi

for i in "${objlst[@]}"; do
    obj_dir="../../asn-ztf/obj_${i}"  # 必须使用双引号
    if cd "${obj_dir}"; then
        echo "Processing source: $i"
        echo "Current directory: $(pwd)"
        python "${baseDirForScriptSelf}/readcsv.py" "$i"
        python "${baseDirForScriptSelf}/calione.py" "$i"
        python "${baseDirForScriptSelf}/plot.py" "$i"
        python "${baseDirForScriptSelf}/wtcali.py" "$i"
        
        # 直接复制文件到目标目录
        if cp -v "${i}.lc" "${lcpath}/" && cp -v "${i}.png" "${lcpath}/"; then
            echo "Files copied successfully to ${lcpath}"
        else
            echo "Error copying files for $i" >&2
        fi
    else
        echo "Error: Cannot enter directory ${obj_dir}" >&2
    fi
done
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "脚本执行时间: $((elapsed_time / 60)) 分钟 $((elapsed_time % 60)) 秒"
echo "脚本执行完毕！"
