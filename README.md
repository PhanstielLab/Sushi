# Sushi

Tools for visualizing genomics data


## Citation

Sushi.r: Flexible, quantitative, and integrative genomic visualizations for publication-quality multi-panel figures.
Phanstiel DH, Boyle AP, Araya CL, Snyder MP. Bioinformatics, In Review. 


## Features

Detailed usage examples are available in the [Vignette](https://github.com/dphansti/Sushi/blob/master/vignettes/Sushi.pdf?raw=true).

## Sushi Installation

<<<<<<< HEAD
1. Install release verson with Bioconductor:
=======
1. Install with Bioconductor:

  Note: R 3.1 is required for installation via Bioconductor. The newest version of R can be downloaded at (www.r-project.org/).  Installation via Bioconductor also requires libcurl and libxml2 which may not be be standard with some Linux distributions and are available for download at (http://curl.haxx.se/libcurl/ and http://xmlsoft.org/, respectively.
>>>>>>> 9093bbff0c29228ed322f0a9b6775763954a7665

  ```
  source("http://bioconductor.org/biocLite.R")
  biocLite("Sushi")
  ```

<<<<<<< HEAD
Note: R 3.1 is required for installation via Bioconductor. The newest version of R can be donwloaded at (www.r-project.org/).  Installation via Bioconductor also requires libcurl and libxml2 which may not be be standard with some Linux distributions and are available for download at (http://curl.haxx.se/libcurl/ and http://xmlsoft.org/, respectviely.


2. Install from source (for previous versions of R >= 2.10):

Download source code from http://www.bioconductor.org/packages/release/bioc/html/Sushi.html

```
R CMD INSTALL Sushi_X.X.X.tar.gz
```

3. The development version of Sushi can be downloaded via Bioconductor

```
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite("Sushi")
```
=======
2. Install for previous versions of R >= 2.10:
  ```
  R CMD INSTALL Sushi_1.0.1.tar.gz
  ```
>>>>>>> 9093bbff0c29228ed322f0a9b6775763954a7665

or via source at http://www.bioconductor.org/packages/devel/bioc/html/Sushi.html


## Contributors

* [Doug Phanstiel](https://github.com/dphansti)
* [Alan Boyle](https://github.com/aboyle)
* Caros Araya

## License
The code is freely available under the GPL (>= 2) license
