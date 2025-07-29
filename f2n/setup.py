#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 使用 setuptools 替代旧的 distutils
from setuptools import setup, find_packages

# 使用 with open(...) 并指定编码，这是更现代和安全的文件读取方式
with open('README.txt', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='f2n',
	# 版本号保持不变
    version='1.2',
	# 使用 find_packages() 自动发现包
    packages=find_packages(),
	# package_data 格式保持不变，确保字体文件被包含
    package_data = {
        'f2n': ['f2n_fonts/*.*']
    },
    # 明确声明运行所需的依赖项
    install_requires=[
        'numpy',
        'astropy',
        'Pillow'
    ],
    description='f2n: A tiny python module to transform FITS images into PNG files.',
    long_description=long_description,
    long_description_content_type='text/plain', # 增加描述类型
    author='Malte Tewes',
    license='GPLv3',
    author_email='malte.tewes[at]epfl.ch',
    url='http://obswww.unige.ch/~tewes/f2n_dot_py/'
)
