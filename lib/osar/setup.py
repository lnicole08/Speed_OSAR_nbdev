from setuptools import setup, find_packages
import os
# Taken from setup.py in seaborn.
# temporarily redirect config directory to prevent matplotlib importing
# testing that for writeable directory which results in sandbox error in
# certain easy_install versions
os.environ["MPLCONFIGDIR"]="."


DESCRIPTION = 'OSAR.'
LONG_DESCRIPTION = """\
Package for analysing the CRITTA output of OSAR experiments.
"""


if __name__ == "__main__":
    setup(
        name='OSAR',
        author='Joses W. Ho',
        author_email='joseshowh@gmail.com',
        maintainer='Joses W. Ho',
        maintainer_email='joseshowh@gmail.com',
        version='0.23.5',
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[
            'numpy>=1.17',
            'scipy>=1.2',
            'pandas>=0.25,!=0.25.2',
            'matplotlib>=3.0',
            'seaborn>=0.9',
            'dabest>=0.3',
            'numba>=0.48'
        ],
        extras_require={'dev': ['pytest>=5.3', 'pytest-mpl>=0.11']},
        python_requires='>=3.5',
        url='http://github.com/ACCLAB/OSAR',
        download_url='http://github.com/ACCLAB/OSAR',
    )
