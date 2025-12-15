---
icon: jet-fighter-up
---

# Quickstart

This page is dedicated to the first steps to start working with SAUNA.

## Installation

Currently, the installation can be done only from source on Linux. If you work with Windows, you can use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) to continue the installation from source.

### Prerequisites

SAUNA relies on the rich Python ecosystem and requires having Python on the system. If the system is fresh, one may install it via:

```bash
sudo apt install python3-dev
```

and make sure to upgrade pip before moving to the [From Source](quickstart.md#from-source) subsection:

```python
pip3 install --upgrade pip
```

Though [NJOY](https://github.com/njoy/NJOY2016) is an exception and has to be installed manually for processing covariance data. One should go to [the official NJOY GitHub repository](https://github.com/njoy/NJOY2016) for further instructions. Still, one shall have git, cmake, and Fortran to properly compile NJOY. They may be installed via the following line:

```bash
sudo apt install git cmake gfortran
```

The other external tool is [AMPX](https://www.ornl.gov/onramp/ampx), developed by Oak Ridge National Laboratory (ORNL), and used here only for reading the COVERX format. Since the usual way to get the data in the format is to get [the SCALE code system](https://scale-manual-does.ornl.gov/), it is assumed that the AMPX is installed alongside SCALE. The other possible way is to get AMPX from [the official ORNL GitLab repository](https://code.ornl.gov/scale/code/scale-public) (though it has not been tested here).

### From Source

Firstly, clone the repository via the following command:

```
git clone https://github.com/AARyzhkov/sauna.git
```

Then, run from the SAUNA root pip install:

```
cd sauna
pip3 install .
```

This concludes the installation process.

An alternative to the `pip install` way is to ensure the package is included in your system's `$PYTHONPATH$` environment variable by adding the appropriate directory. The package can also be added directly in your Python code:

```python
import sys
sys.path.append('./path/to/the/package/sauna')
```

In addition, one has to install the following packages manually if `pip install` is not used:

* [SANDY](https://github.com/luca-fiorito-11/sandy)
* [serpentTools](https://github.com/CORE-GATECH-GROUP/serpent-tools)
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
