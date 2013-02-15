from distutils.core import setup

import SARlabpy
README = os.path.join(os.path.dirname(__file__), 'README.txt')
long_description = open(README).read() + 'nn'

setup(name='SARlabpy',
      version=SARlabpy.version,
      author='Stefan A Reinsberg and other SARlab members',
      url='http://pfeifer.phas.ubc.ca/SARlabpy',
      download_url='http://pfeifer.phas.ubc.ca/SARlabpy/files',
      description='Data Analysis of MR data acquired on BRUKER and other systems',
      long_description=long_description, 
      Iackage_dir={'': 'SARlabpy'},
      py_modules=['SARlabpy'],
      provides=['SARlabpy'],
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
