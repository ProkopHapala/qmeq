
import os
import sys
import numpy as np

from setuptools import setup

def get_ext_modules():
    return []

long_description = open('README.rst').read()

setup(
      packages=['qmeq',
                'qmeq/approach',
                'qmeq/approach/base',
                'qmeq/approach/elph',
                'qmeq/builder',
                'qmeq/specfunc',
                'qmeq/tests',
                'qmeq/wrappers',],
      package_data={}, # Removed package_data to avoid including .pyx and .c files
      zip_safe=False,
      include_dirs=[np.get_include()],
      ext_modules=get_ext_modules())
