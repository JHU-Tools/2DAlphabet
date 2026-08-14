#!/usr/bin/env python

import setuptools, sys

requires = [
      'pandas',
      'Pillow',
      'tabulate'
]
if sys.version_info.major == 3:
    requires.extend(['doxypypy','pytest','pytest-cov'])
    
setuptools.setup(name='TwoDAlphabet',
      version='2.0',
      description='Framework for performing 2D binned-likelihood fits with one background derived from a transfer function.',
      author='Matej Roguljic',
      author_email='matej.roguljic@gmail.com',
      url='https://gitlab.cern.ch/jhu-tools/2DAlphabet',
      packages=setuptools.find_packages(),
      include_package_data=True,
      # cmdclass={'install': AddToPath},
      install_requires = requires
)