# Sushi

Tools for visualizing genomics data

## Citation

Phanstiel DH, Boyle AP, Araya CL, Snyder MP. Sushi.R: flexible, quantitative
and integrative genomic visualizations for publication-quality multi-panel
figures. Bioinformatics. 2014 Oct;30(19):2808-10.

## Features

Detailed usage examples are available in the [Vignette](https://github.com/dphansti/Sushi/blob/master/vignettes/Sushi.pdf?raw=true).

## Sushi Installation

1. To install the latest version directly from github open R and run the following commands (requires the R package devtools which is available through CRAN):

```
library("devtools")
install_github("dphansti/Sushi")
```

2. Install release verson with Bioconductor:

 ```
 source("http://bioconductor.org/biocLite.R")
 biocLite("Sushi")
 ```

 Note: R 3.1 is required for installation via Bioconductor. The newest version of R can be downloaded at (www.r-project.org/).  Installation via Bioconductor also requires libcurl and libxml2 which may not be be standard with some Linux distributions and are available for download at (http://curl.haxx.se/libcurl/ and http://xmlsoft.org/, respectively.


3. Install from source (for previous versions of R >= 2.10):

 Download source code from http://www.bioconductor.org/packages/release/bioc/html/Sushi.html

 ```
 R CMD INSTALL Sushi_X.X.X.tar.gz
 ```

4. The development version of Sushi can be downloaded via Bioconductor

 ```
 source("http://bioconductor.org/biocLite.R")
 useDevel()
 biocLite("Sushi")
 ```

 or via source at http://www.bioconductor.org/packages/devel/bioc/html/Sushi.html


## Contributors

* [Doug Phanstiel](https://github.com/dphansti)
* [Alan Boyle](https://github.com/aboyle)
* Caros Araya

## License
The code is freely available under the GPL (>= 2) license
