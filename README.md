# cohomCalg
[![GitHub release](https://img.shields.io/github/release/BenjaminJurke/cohomCalg.svg)](https://github.com/BenjaminJurke/cohomCalg/releases) [![GitHub issues](https://img.shields.io/github/issues/BenjaminJurke/cohomCalg.svg)](https://github.com/BenjaminJurke/cohomCalg/issues) [![GitHub license](https://img.shields.io/badge/license-GPL%203-lightgrey.svg)](https://raw.githubusercontent.com/BenjaminJurke/cohomCalg/master/LICENSE)

> ![cohomCalg logo with Koszul extension](https://raw.githubusercontent.com/BenjaminJurke/cohomCalg/master/cohomCalg_koszul_logo.png)  
> A software package for computation of sheaf cohomologies  
> for line bundles on toric varieties.

**Authors:** [Ralph Blumenhagen](http://wwwth.mpp.mpg.de/members/blumenha/), [Benjamin Jurke](https://benjaminjurke.net), [Thorsten Rahn](http://thorsten-rahn.net), Helmut Roschy

# Project history

The algorithm for the computation of sheaf cohomologies for line bundles on toric varieties presented in [arXiv:1003.5217 [hep-th]](http://arxiv.org/abs/1003.5217) *"Cohomology of Line Bundles: A Computational Algorithm"* has been implemented in a convenient and high-performance C/C++ application called **cohomCalg**. The optional **cohomCalg Koszul** extension serves as a Mathematica 7 frontend and allows for the easy computation of hypersurface and complete intersection cohomologies, following the material presented in [arXiv:1010.3717 [hep-th]](http://arxiv.org/abs/1010.3717).

About three months after the initial conjecture's preprint release, a proof of the algorithm was presented in [arXiv:1006.2392 [hep-th]](http://arxiv.org/abs/1006.2392), which also clarifies much of the underlying mathematical structures. At the same time an independant proof was developed in [arXiv:1006.0780 [math.AG]](http://arxiv.org/abs/1006.0780) - published in fact a few days earlier - which utilized alternative methods. 

The [full background story](https://cohomcalg.benjaminjurke.net/) can be found here along with some simple examples.

The source code is freely available under the GNU GPL v3 license terms. Furthermore, the implementation makes use of The Polyhedral Library (or [PolyLib](http://icps.u-strasbg.fr/polylib/) for short), which is available here under the same license.

# Documentation

A full documentation and several example input files are included in the cohomCalg download package, see "[manual.pdf](https://github.com/BenjaminJurke/cohomCalg/blob/master/manual.pdf)". Furthermore, you can [eMail us](mailto:mail@benjaminjurke.net?subject=cohomCalg) for technical support or other related questions. For ease of use the package includes pre-compiled binaries for Microsoft Windows both in 32- and 64-bit versions, with the latter being considerably faster.


# Changelog

* **v0.31c** (January 14, 2017): Minor update due to GCC 6.2 compatibility issues, recompiled Windows binaries.
* **v0.31b** (April 18, 2012): Minor bugfix due to GCC 4.6 compatibility issues.
* **v0.31** (May 25, 2011): Multi-core support, new vector bundle routines, overall improvements.
* **v0.21** (October 18, 2010): cohomCalg Koszul extension added!
* **v0.13** (July 23, 2010): Integration mode added, see command line option `--integrated` in manual.
* **v0.12** (June 25, 2010): Several minor bugfixes and improvements.
* **v0.11** (May 4, 2010): Original public release of the cohomCalg implementation.
* **v0.04** (March 29, 2010): Original public release of the Mathematica 7 script, which is still available [here](https://github.com/BenjaminJurke/cohomCalg/tree/master/old/cohomcalg-script-v004).


# Known bugs & Shortcomings

Certain invalid input data is at the moment not safely handled by the PolyLib. For the moment, the crash situation is circumvented in a Quick&Dirty manner. In such cases you will see a message like `Counting of the rationoms errorneous - is your input geometry valid?` in those cases, which caused the original release version to crash. So far, all such cases could be traced back to bad input data. Further implementation-related shortcomings are explained in the manual.


# License & Proper Citation

As mentioned, the entire package is published under the GNU GPL v3 License, as required by the included PolyLib. This means that any derivative work also has to be published under the GPL v3 License or an equivalent license. The package contains a small text file "Proper Citation.txt" providing a BibTeX entry for cohomCalg, which you can use if you use the program in your work. 

Or you can simply **Copy&Paste the following BibTeX snippet** provided for convenience:


```tex
@Article{Blumenhagen:2010pv,
   author    = "Blumenhagen, Ralph and Jurke, Benjamin 
                and Rahn, Thorsten and Roschy, Helmut",
   title     = "{Cohomology of Line Bundles: A Computational Algorithm}",
   journal   = "J. Math. Phys.",
   volume    = "51",
   pages     = "103525",
   issue     = "10",
   year      = "2010",
   doi       = "10.1063/1.3501132",
   eprint    = "1003.5217",
   archivePrefix = "arXiv",
   primaryClass  = "hep-th"}

@Misc{cohomCalg:Implementation,
   title     = "{cohomCalg package}",
   howpublished  = "Download link",
   url       = "http://wwwth.mppmu.mpg.de/members/blumenha/cohomcalg/",
   note      = "High-performance line bundle cohomology computation based on \cite{Blumenhagen:2010pv}",
   year      = "2010"}
````


# Related Links

In order to derive the Stanley-Reisner ideal, which is a required input for the program, you may want to take a look at [TOPCOM](http://www.rambau.wm.uni-bayreuth.de/TOPCOM/), which can also enumerate all possible fans for a given set of vertices. The Maple script package [SCHUBERT](http://stromme.uib.no/schubert/) can be used to compute intersection numbers and further geometrical quantities of toric varieties. Furthermore, there is the package [PALP](http://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html) which is useful for computing invariants of hypersurfaces, Mori cone vectors etc. You may also want to take a look at the [SAGE Library](http://www.sagemath.org/) of freely available mathematical software. The [Macaulay2](http://www.math.uiuc.edu/Macaulay2/) software, which allows similar computations, was heavily used during the development process.