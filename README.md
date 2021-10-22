# Wphase 

The Wphase package is an extensive source-inversion tool allowing rapid characterization of seismic sources and is a reliable and straightforward method to resolve first-order attributes of large earthquakes.   

For a complete documentation, please refer to [the Wphase wiki page](http://wphase.unistra.fr/wiki/doku.php/wphase)

Main contributors are Z. Duputel, L. Rivera and H. Kanamori (California Institute of Technology and Universite de Strasbourg / CNRS).

Wphase is licensed under the GNU General Public License (GPL) v3.


[![GPL](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl.html) 
[![Build Status](https://travis-ci.org/eost/wphase.svg?branch=master)](https://travis-ci.org/eost/wphase)

## Citation

If you use the W-phase package for your own research, please cite the following articles:

* Z. Duputel, L. Rivera, H. Kanamori, and G. Hayes, 2012. W-phase fast source inversion for moderate to large earhquakes (1990 - 2010), Geophysical Journal International, v. 189, iss. 2, p. 1125-1147.
* G.P. Hayes, L. Rivera   and H. Kanamori, 2009. Source inversion of the W phase: real-time implementation and extension to low magnitudes, Seismol. Res. Let., v. 3, 800-805.
* H. Kanamori and L. Rivera, 2008. Source inversion of W phase: speeding tsunami warning, Geophys. J. Int., v. 175, 222-238.

## Installation

### Dependencies
* csh shell
* gcc, gfortran
* python2 or python3
* You have to install numpy, matplotlib, basemap and netCDF4 to run some python scripts which make figures.
* rdseed

To install these dependencies on MacOs, [you can refer to this page](http://wphase.unistra.fr/wiki/doku.php/wphase:macos).

### Building the code

To install the code, we must first setup a few environment variables. If you use csh or tcsh:

```
setenv RDSEED       /path/to/rdseed/executable
setenv GF_PATH      /path/to/greens/functions/database
setenv WPHASE_HOME  /path/to/wphase/package
```

If you use bash:

```
export RDSEED=/path/to/rdseed/executable
export GF_PATH=/path/to/greens/functions/database
export WPHASE_HOME=/path/to/wphase/package
```

These variables are necessary both at installation time and at run time. In addition to these variables, it is handy to include the wphase bin directory into the PATH environment variable:

```
setenv PATH         ${PATH}:$WPHASE_HOME/bin
```

or in bash:

```
export PATH=${PATH}:$WPHASE_HOME/bin
```

All theses variable assignment can be included in your .tcshrc or your .bashrc file.

Once you have downloaded the last version of the W-phase package and that the above environment variables are properly defined. Then you can proceed as follows to compile the package:

```
cd ${WPHASE_HOME}/src
make -B
```

## How to run W-phase

Instructions on how to run and use Wphase are available at the [Wphase wiki page](http://wphase.unistra.fr/wiki/doku.php/wphase). 

There is a quite complete document in the [Wphase documentation](http://wphase.unistra.fr/wiki/doku.php/wphase:documentation) and you can also check out the [Wphase tutorial page](http://wphase.unistra.fr/wiki/doku.php/wphase:tutorial) which includes an example of run.

## Development
Development is hosten on GitHub in the [eost/wphase repository](https://github.com/eost/wphase).

-- 

**Report bugs to: <zacharie.duputel@unistra.fr> and <luis.rivera@unistra.fr>**
