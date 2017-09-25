---
author: Richard Glennie
generator: pandoc
title: 'HMM Madness: R, ARMA, and TMB'
viewport: 'width=device-width, initial-scale=1'
---

<div class="container-fluid main-container">

<div id="header" class="fluid-row">

HMM Madness: R, ARMA, and TMB {#hmm-madness-r-arma-and-tmb .title .toc-ignore}
=============================

#### *Richard Glennie* {#richard-glennie .author}

#### *21/09/2017* {#section .date}

</div>

Code to fit Poisson Hidden Markov Model (see REF) in R, RcppAramdillo
and Template Model Builder. Coding of the algorithm is strikingly
similar between the three programming approaches. Speed comparisons
showed TMB to be faster in all cases, particularly faster when sample
size or the number of hidden states was high. Motivation for this
comparison was to consider best approach for obtaining non-parametric
bootstrap estimates of variance, requiring the HMM to be fit many times.

A basic knowledge of C++ makes coding the HMM algorithm in RcppArmadillo
just as easy as coding in R. Implementing Template Model Builder
required more knowledge of C++ and was less easy to use (see some notes
on TMB below).

<div id="file-descriptions" class="section level3">

### File descriptions

-   `poisson_hmm.R`: R functions to simulate Poisson HMM data, transform
    between natural and working parameters, and functions to fit HMM
    using R, RcppArmadillo, or Template Model Builder.
-   `hmm_arma.cpp`: C++ code to compute HMM likelihood using Armadillo
    and Rcpp.
-   `hmm_tmb.cpp`: C++ code to compute HMM likelihood using Template
    Model Builder.
-   `simulation.R`: R function for performing simulation study on
    Poisson HMM using three programming approaches.
-   `run.R`: R code to run a simulation study and output the average
    time taken to fit 999 data sets with R, ARMA and TMB under a few
    different sample sizes and hidden state dimensions.

</div>

<div id="speed-comparison" class="section level3">

### Speed Comparison

  n       NStates   R       ARMA    TMB
  ------- --------- ------- ------- ------
  100     2         0.04    0.00    0.01
  1000    2         0.25    0.03    0.02
  10000   2         2.73    0.29    0.18
  100     4         0.80    0.10    0.02
  1000    4         6.14    1.02    0.14
  10000   4         62.48   10.51   1.61
  100     6         2.33    0.43    0.03
  1000    6         17.80   4.25    0.27
  10000   6         NA      51.49   3.37

  : Average time (in secs) taken to maximise Poisson HMM likelihood,
  with different numbers of states (NStates), for 99 data sets, each
  with sample size n, using R, RcppArmadillo (ARMA) and Template Model
  Builder (TMB) (NA means &gt; 5 minutes)

</div>

<div id="notes-on-using-tmb" class="section level3">

### Notes on Using TMB

-   A familiarity with C++ templates isn’t essential, but is useful for
    debugging. The automatic differentation relies on the correct
    variables being declared as `<Type>`.
-   The TMB library is not as robust as Eigen; sometimes operations that
    work in Eigen do not work with TMB `matrix` and `vector` objects,
    despite being based on Eigen. In most cases, this is because no
    polymorphic alternative has been written to handle CppAD types. For
    example, `(v % w) * A` where `v`,`w` are row vectors and `A` is a
    matrix, fails because TMB does not have the required support for
    operation chains.
-   The dynamic library produced when compiling TMB objects has to be
    loaded into R. On some platforms (including mine), if you re-compile
    the object and re-load the dynamic library, R either returns an
    error (because conflicting libraries are loaded) or simply ignores
    the re-load and retains the old version of the dynamic library. This
    is particularly frustrating when you are developing code: edit TMB
    function, re-compile, re-load, does not re-load, cannot run editted
    function and debug. The only solution I found was to restart R every
    time I wanted to re-load the dynamic library.

</div>

</div>
