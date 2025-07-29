# -*- coding: utf-8 -*-
import os
import subprocess
import tempfile

# 图像文件夹和坐标文件
image_dir = 'alipy_out'  # 图像文件夹路径
pos_file = 'PHL1092.pos'  # 你的星点坐标文件
log_dir = 'psf_logs'  # 日志输出文件夹

# 创建日志文件夹
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

def run_psfmeasure_subprocess(image_path, pos_file, log_path):
    """使用subprocess直接调用IRAF命令"""
    
    # 创建临时的IRAF脚本
    script_content = f"""noao
obsutil
set stdimage=imt1024
set imcur=""
set imtype=fits
psfmeasure {image_path} imagecur={pos_file} coords=markall size=FWHM logfile={log_path}
logout
"""
    
    # 写入临时脚本文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cl', delete=False) as f:
        f.write(script_content)
        script_file = f.name
    
    try:
        # 设置环境变量
        env = os.environ.copy()
        env['IRAF_IMCUR'] = ''
        env['DISPLAY'] = ':0.0'
        env['TERM'] = 'xterm'
        
        # 运行IRAF脚本
        result = subprocess.run(['cl', '<', script_file], 
                              env=env, 
                              capture_output=True, 
                              text=True,
                              timeout=300)  # 5分钟超时
        
        if result.returncode != 0:
            print(f"Error processing {image_path}: {result.stderr}")
        else:
            print(f"Successfully processed {image_path}")
            
    except subprocess.TimeoutExpired:
        print(f"Timeout processing {image_path}")
    except Exception as e:
        print(f"Exception processing {image_path}: {e}")
    finally:
        # 清理临时文件
        if os.path.exists(script_file):
            os.unlink(script_file)

def main():
    for fname in os.listdir(image_dir):
        if fname.endswith('.fits'):
            image_path = os.path.join(image_dir, fname)
            log_path = os.path.join(log_dir, fname.replace('.fits', '_psf.txt'))
            print(f"Processing {fname}...")
            run_psfmeasure_subprocess(image_path, pos_file, log_path)
    print("PSF measurement completed for all images.")

if __name__ == '__main__':
    main()
