Code to fit Poisson Hidden Markov Model (see Zucchini et al. 2016) in R, RcppAramdillo, Template Model Builder, 
and Python. 

Coding of the algorithm is strikingly
similar between the all programming approaches. 

A speed comparison was performed between R, RcppArmadillo and Template Model Builder. 
TMB was faster in all cases, particularly faster when sample
size or the number of hidden states was high. Motivation for this
comparison was to consider best approach for obtaining non-parametric
bootstrap estimates of variance, requiring the HMM to be fit many times.
A basic knowledge of C++ makes coding the HMM algorithm in RcppArmadillo
just as easy as coding in R. Implementing Template Model Builder
required more knowledge of C++ and was less easy to use (see some notes
on TMB below).


### File descriptions

-   <code>poisson\_hmm.R</code>: R functions to simulate Poisson HMM
    data, transform between natural and working parameters, and
    functions to fit HMM using R, RcppArmadillo, or Template Model
    Builder.
-   <code>hmm\_arma.cpp</code>: C++ code to compute HMM likelihood using
    Armadillo and Rcpp.
-   <code>hmm\_tmb.cpp</code>: C++ code to compute HMM likelihood using
    Template Model Builder.
-   <code>py_hmm.py</code>: python3 code to fit poisson HMM 
-   <code>simulation.R</code>: R function for performing simulation
    study on Poisson HMM using three programming approaches.
-   <code>run.R</code>: R code to run a simulation study and output the
    average time taken to fit 999 data sets with R, ARMA and TMB under a
    few different sample sizes and hidden state dimensions.

### Speed Comparison

<table>
<caption>Average time (in secs) taken to maximise Poisson HMM likelihood, with different numbers of states (NStates), for 99 data sets, each with sample size n, using R, RcppArmadillo (ARMA) and Template Model Builder (TMB) (NA means &gt; 5 minutes)</caption>
<thead>
<tr class="header">
<th align="right">n</th>
<th align="right">NStates</th>
<th align="right">R</th>
<th align="right">ARMA</th>
<th align="right">TMB</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">100</td>
<td align="right">2</td>
<td align="right">0.04</td>
<td align="right">0.00</td>
<td align="right">0.01</td>
</tr>
<tr class="even">
<td align="right">1000</td>
<td align="right">2</td>
<td align="right">0.25</td>
<td align="right">0.03</td>
<td align="right">0.02</td>
</tr>
<tr class="odd">
<td align="right">10000</td>
<td align="right">2</td>
<td align="right">2.73</td>
<td align="right">0.29</td>
<td align="right">0.18</td>
</tr>
<tr class="even">
<td align="right">100</td>
<td align="right">4</td>
<td align="right">0.80</td>
<td align="right">0.10</td>
<td align="right">0.02</td>
</tr>
<tr class="odd">
<td align="right">1000</td>
<td align="right">4</td>
<td align="right">6.14</td>
<td align="right">1.02</td>
<td align="right">0.14</td>
</tr>
<tr class="even">
<td align="right">10000</td>
<td align="right">4</td>
<td align="right">62.48</td>
<td align="right">10.51</td>
<td align="right">1.61</td>
</tr>
<tr class="odd">
<td align="right">100</td>
<td align="right">6</td>
<td align="right">2.33</td>
<td align="right">0.43</td>
<td align="right">0.03</td>
</tr>
<tr class="even">
<td align="right">1000</td>
<td align="right">6</td>
<td align="right">17.80</td>
<td align="right">4.25</td>
<td align="right">0.27</td>
</tr>
<tr class="odd">
<td align="right">10000</td>
<td align="right">6</td>
<td align="right">NA</td>
<td align="right">51.49</td>
<td align="right">3.37</td>
</tr>
</tbody>
</table>

### Notes on Using TMB

-   A familiarity with C++ templates isn't essential, but is useful for
    debugging. The automatic differentation relies on the correct
    variables being declared as <code>&lt;Type&gt;</code>.
-   The TMB library is not as robust as Eigen; sometimes operations that
    work in Eigen do not work with TMB <code>matrix</code> and
    <code>vector</code> objects, despite being based on Eigen. In most
    cases, this is because no polymorphic alternative has been written
    to handle CppAD types. For example, <code>(v % w) \* A</code> where
    <code>v</code>,<code>w</code> are row vectors and <code>A</code> is
    a matrix, fails because TMB does not have the required support for
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
