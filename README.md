# SAUNA

[![Documented on GitBook](https://img.shields.io/static/v1?message=Documented%20on%20GitBook&logo=gitbook&logoColor=ffffff&label=%20&labelColor=5c5c5c&color=3F89A1)](https://sauna-dev.gitbook.io/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/AARyzhkov/sauna/blob/master/LICENSE)

A Python package for sensitivity and uncertainty analysis

## Description

This package has been developed for conducting seamless systematic sensitivity and uncertainty analysis with respect to nuclear data. It provides capabilities from producing covariance nuclear data to propagating the uncertainties and finding target accuracy requirements (TARs) for nuclear data to satisfy TARs for reactor functionals (responses). SAUNA allows one to systematically analyze nuclear data influence on functionals leveraging Python capabilities.

Currently, the documentation is hosted via [gitbook](https://sauna-dev.gitbook.io).

## Installation

Currently, the installation can be done only from source on Linux. If you work with Windows, you can use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) to continue the installation from source.

### Prerequisites

SAUNA relies on the rich Python ecosystem and requires having Python on the system. If the system is fresh, one may install it via:

```bash
sudo apt install python3-dev
```

and make sure to upgrade pip before moving to the [From Source](README.md#from-source) subsection:

```python
pip3 install --upgrade pip
```

### External tools

Though [NJOY](https://github.com/njoy/NJOY2016) is an exception and has to be installed manually for processing covariance data. One should go to [the official NJOY GitHub repository](https://github.com/njoy/NJOY2016) for further instructions. Still, one shall have git, cmake, and Fortran to properly compile NJOY. They may be installed via the following line:

```bash
sudo apt install git cmake gfortran
```

After NJOY is installed, make sure a path to the NJOY executable is provided in the following manner:

```shellscript
export NJOY='/home/your-username/NJOY2016/build/njoy'
```

The other external tool is [AMPX](https://www.ornl.gov/onramp/ampx), developed by Oak Ridge National Laboratory (ORNL), and used here only for reading the COVERX format. Since the usual way to get the data in the format is to get [the SCALE code system](https://scale-manual-does.ornl.gov/), it is assumed that the AMPX is installed alongside SCALE. The other possible way is to get AMPX from [the official ORNL GitLab repository](https://code.ornl.gov/scale/code/scale-public) (though it has not been tested here).

### From Source

Firstly, clone the repository via the following command:

```
git clone https://github.com/AARyzhkov/sauna.git
```

#### via pip

Then, run from the SAUNA root pip install (there is a point at the end of the second line):

```
cd sauna
pip3 install .
```

This concludes the installation process.

#### via path

An alternative to the ```pip install``` way is to ensure the package is included in your system's ```$PYTHONPATH$``` environment variable by adding the appropriate directory. The package can also be added directly in your Python code:

```python
import sys
sys.path.append('./path/to/the/package/sauna')
```

In addition, one has to install the following packages manually if ```pip install``` is not used:

* [SANDY](https://github.com/luca-fiorito-11/sandy) for providing the interface for NJOY
* [serpentTools](https://github.com/CORE-GATECH-GROUP/serpent-tools) for providing API for interacting with Serpent results
* [NumPy](https://github.com/numpy/numpy)
* [SciPy](https://github.com/scipy/scipy)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [pandas](https://github.com/pandas-dev/pandas)

This can be done in a batch via the following line:
```python
pip3 install sandy serpentTools numpy scipy matplotlib pandas
```

One also may install texlive to make it possible to use TeX for rendering:
```bash
sudo apt install texlive texlive-latex-extra texlive-fonts-recommended dvipng cm-super
```

## License

SAUNA is distributed under the MIT/X license.
