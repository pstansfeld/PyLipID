# setup.py
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.rst').read_text(encoding='utf-8')

# read version info
import re
VERSIONFILE="pylipid/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

# setup
setup(
    name='pylipid',
    version=verstr,
    description='PyLipID - A Python Library For Lipid Interaction Analysis',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/pstansfeld/PyLipID',
    author='Wanling Song',
    author_email='wanling.song@hotmail.com',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    keywords='simulation tools, network community, binding site',
    python_requires='>=3.9',
    packages=find_packages(),
    install_requires=[
        "mdtraj>=1.9.7",
        "numpy>=1.21",
        "pandas>=1.3",
        "matplotlib>=3.5",
        "networkx>=2.6",
        "scipy>=1.7",
        "logomaker>=0.8",
        "statsmodels>=0.13",
        "scikit-learn>=1.0",
        "tqdm>=4.50",
        "p_tqdm>=1.3",
        "kneebow>=0.1",
    ]
)
