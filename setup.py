from setuptools import setup, find_packages
from os import path
here = path.abspath(path.dirname(__file__))

def get_version():
    with open(path.join(here, "nanomonsnv/version.py")) as hin:
        for line in hin:
            if line.startswith("__version__"):
                version = line.partition('=')[2]
                return version.strip().strip('\'"')
    raise ValueError('Could not find version.')

setup(
      name='nanomonsnv',
      version=get_version(),
      description="Python programs for Detecting somatic SNV in nanopore sequensing data.",
      long_description="""""",

      classifiers=[
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 3 - Alpha',
          # Indicate who your project is intended for
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      
      keywords='Bio-informatics',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gmail.com',
      url='https://github.com/ncc-ccat-gap/nanomonsnv.git',
      license='GPL3',
      
      packages = find_packages(exclude = ['tests']),
      install_requires=[
      ],
      entry_points = {'console_scripts': ['nanomonsnv = nanomonsnv:main']},
      test_suite = 'unit_tests.suite'
)
