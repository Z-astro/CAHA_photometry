#!/usr/bin/env bash

# 重要！必须加载 conda 初始化
source ~/miniconda3/etc/profile.d/conda.sh
conda activate astroconda
echo 'activate conda env: astroconda (python==3.7.16)'


curpath=`pwd`
# rootdir=`dirname ${curpath}`
# rootdir=/media/yao/cahayao2
rootdir=/mnt/WDsata/CAHA
# baseDirForScriptSelf=$(cd "$(dirname "$0")"; pwd)
baseDirForScriptSelf=/mnt/WDsata/CAHA/codes/step1

cd ${rootdir}
echo '***cd '${rootdir}'/***'
echo 
if ! [ -d objects ];then
  mkdir objects
  echo '***mkdir objects/***'
  echo 
fi

objdir=/mnt/WDsata/CAHA/objects/
# objdir=/media/yao/cahayao2/object_workyao/object/
# echo ***${objdir}***
cd ${curpath}
echo '***cd '${curpath}'/***'
echo ''

update_reduce=`cat update_reduce.lst`
# echo ${update_reduce}
for dirnm in `ls *CAFOS -d`;do 
  result=$(echo $update_reduce | grep "${dirnm}")
  if [[ "${result}" == "" ]];then

    echo '***start check_fits.py***'
    # python8 ${baseDirForScriptSelf}/check_fits.py ${dirnm}
    python ${baseDirForScriptSelf}/check_fits.py ${dirnm}
    echo '***end   check_fits.py***'
    echo ''

    cd ${curpath}/${dirnm}
    echo '***cd '${dirnm}'/***'
    echo ''

    date=${dirnm:0:6}
    for photodir in `ls ${date}_photo* -d`;do
      rm -rf ${photodir}
      echo '***rm -rf '${photodir}'***'
      echo ''
    done

    echo '***start getobs_cp2phot.py***'
    # python8 ${baseDirForScriptSelf}/getobs_cp2phot.py
    python ${baseDirForScriptSelf}/getobs_cp2phot.py
    echo '***end   getobs_cp2phot.py***'
    echo ''

    for photodir in `ls *_photo* -d`;do
      cd ${photodir}
      echo '***cd '${photodir}'/***'
      echo ''

      echo '***start reduce_imgpy3.py***'
      python ${baseDirForScriptSelf}/reduce_imgpy3.py
      echo '***end   reduce_imgpy3.py***'
      echo ''

      echo '***start getobjlst_cp2img.py***'
      # python8 ${baseDirForScriptSelf}/getobjlst_cp2img.py ${objdir}
      python ${baseDirForScriptSelf}/getobjlst_cp2img.py ${objdir}
      echo '***end   getobjlst_cp2img.py***'
      echo ''

      cd ..
      echo '***cd ../***'
      echo ''
    done
    cd ..
    echo '***cd ../***'
    echo ''
    echo ${dirnm} >> update_reduce.lst
  else
    echo '***'${dirnm}' reduced***'
  fi
done
