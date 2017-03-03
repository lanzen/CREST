from setuptools import setup, find_packages
import os

version = '2.0'

setup(name='LCAClassifier',
      version=version,
      description="LCAClassifier",
      long_description=open("README.txt").read(),
      # Get more strings from
      # http://pypi.python.org/pypi?:action=list_classifiers
      classifiers=[
        "Programming Language :: Python",
        ],
      keywords='',
      author='',
      author_email='',
      url='https://github.com/lanzen/CREST',
      license='GPL',
      packages=find_packages('src', exclude=['ez_setup']),
      package_dir={'': 'src'},
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
          'numpy',
          'biopython',
          'biom-format',
          ],
      entry_points={
                    'console_scripts': [
                              'classify = LCAClassifier.classify:main',
                              ]
                    },
      )
