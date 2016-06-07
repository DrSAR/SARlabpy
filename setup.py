from setuptools import setup
import os

#import sarpy
from sarpy import version
from sarpy import release
README = os.path.join(os.path.dirname(__file__), 'README.rst')
long_description = open(README).read() + 'nn'

# insidious problem here:
# if this is not wrapped in a if __name__ ... then the sphinx harvesting script
# barfs with an obstruse error
# option -b not recognized
if __name__ == '__main__':
    setup(name='sarpy',
      version=version,
      author='Stefan A Reinsberg and other SARlab members',
      url='https://github.com/DrSAR/SARlabpy',
      description='Data Analysis of MR data acquired on BRUKER and other systems',
      long_description=long_description, 
      package_dir={'': 'sarpy'},
      py_modules=['sarpy','BRUKERIO'],
      provides=['sarpy'],
      keywords='MRI image processing data analysis',
      license='General Public License v3',
      classifiers=['Development Status :: 2 - Pre-Alpha',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2',
                   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                   'Topic :: Image Processing',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.',
                  ],
     )
