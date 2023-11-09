---
output:
  html_document: default
  pdf_document: default
---
# Changes and Updates for PracTools package

# PracTools 1.4.1

*    Incorrect calculation of 1-draw selection probabilities corrected in vignette, Distance-and-MOS-PSUs.

*    Added check that 1-draw selection probabilities in parameter pp sum to 1 in the functions, BW2stagePPS, BW3stagePPs, deff, and deffS.

*    Added a new vignette titled "Design Effects and Effective Sample Size".

*    deffCR output updated to include stratum estimates of rho and coefficients of variation of weights

# PracTools 1.4

*    New functions,  CompMOS and GeoMinMOS, added. Edited clusOpt2 to allow vector input for some parameters. Vignette "Selection of Appropriate PracTools Sample Size Function" added.
    
*    CompMOS computes a composite measure of size variable for domain-based sampling that accounts for desired sampling rates and anticipated response rates of domain units.
  
*    GeoMinMOS takes geographic areas as input, determines whether any are less than a specified minimum measure of size, and combines areas that are too small with other nearby areas. 
    
*    There are 19 different functions in PracTools for computing sample sizes. The vignette "Selection of Appropriate PracTools Sample Size Function" gives some guidance on how to select a function for different applications.

# PracTools 1.3

*    New functions added are: GeoDistPSU, GeoDistMOS, and nContOpt. BW2stagePPS and BW2stagePPSe updated to handle minor issues with certainty selections.

*    BW2stagePPS and BW2stagePPSe updated to deal with PSUs whose input selection probability is 1. In BW2stagePPS such PSUs are excluded from calculations, results are produced based on PSUs with pp < 1, and a warning is printed with the number of PSUs that were excluded.  In BW2stagePPSe PSUs with pp = 1 are not allowed. If the input contains any PSUs with pp >= 1, the function is halted.  These PSUs should be removed before calling BW2stagePPSe.

*    GeoDistPSU was added that combines geographic units into PSUs based on latitude and longitude of their centroids.  A maximum within-PSU distance can be set to control the size of PSUs.

*    GeoDistMOS is a companion function to GeoDistPSU. Geographic units are formed into PSUs based on latitude and longitude of their centroids as in GeoDistPSU.  A maximum within-PSU distance can be set to control the size of PSUs.  In addition, a measure of size (MOS) for each PSU is used to limit the maximum MOS allowed for a PSU.

*    A vignette was also added that describes how to use GeoDistPSU and GeoDistMOS. When the geographic units are in the US, the plot_usmap function in ggplot2 can be used to display a map of the PSUs that were formed.  An html version of the vignette can be viewed with
```
   vignette("Distance-and-MOS-PSUs", package="PracTools")
```

*    nContOpt computes the sample size required to estimate the mean of a continuous variable by optimizing the numbers of take-alls (certainties) and non-take-all units selected by probability sampling.  The sample size is based on splitting the sample between take-alls and non-take-alls in a way that achieves either a target coefficient of variation or a target variance for an estimated mean. The sample design for the non-take-alls can be either simple random sampling or probability proportional to size sampling.

# PracTools 1.2.8

*    nAuditMUS function added for monetary unit sampling in audit applications. nAuditAttr was edited to produce somewhat larger sample sizes.  Users will notice a difference from the previous version mainly when the population sizes are small. 

# PracTools 1.2.7

*    Edited help file for NRFUopt to better describe contents of the output; removed “Expected total cases (2‐phase)” from output. Edited clusOpt2 to allow only scalar parameter inputs. 

# PracTools 1.3
*    Corrected error in BW3stagePPSe. (browser() inadvertently left in function.) 
    
# PracTools 1.2.5

*    Corrected error in deffCR in estimation of $\small {{\rho }_{h}}$ for stratified samples. 
 
# PracTools 1.2.4

*    The function nAuditAttr was added to calculate sample sizes for attribute samples in an audit of financial records. 
   
# PracTools 1.2.3

*    Modified BW2stagePPS, BW2stagePPSe, BW2stageSRS, BW3stagePPS, and BW3stagePPSe to not return NA for the SSU variance component when a PSU contains only 1 SSU.   

*    Parameters lonely.SSU and lonely.TSU were added. With the default values of “mean”, any NAs are replaced with the mean of non-missing values. The other allowable value is “zero” in which case any missing value is replaced by 0. 

*    In help files, updated references from sections in VDK (2013) to VDK (2018). 
 
# PracTools 1.2.2

*    Corrected error in v.1.2.1 warning in nContMoe.  The warning should be for the combination: moe.sw=2, e >= 1. 

# PracTools 1.2.1

*    nContMoe edited to eliminate restriction that parameter e (margin of error) must be less than 1.  If moe.sw = 2 and e >= 1, a warning is printed:  
```    
   WARNING: e >= 1. This parameter setting leads to the lower limit of a 2-sided alpha level confidence interval being 0 or less. You may want to adjust the value of e. 
```

# PracTools 1.1

*    A population called ThirdGrade of 2,427 third graders’ math and science scores was added. The data are from the Third International Mathematics and Science Study.  
    
*    In some applications the NRadjClass function would fail if there were many duplicate values of predicted probabilities generated by the pclass function.  pclass was modified to jitter the predictions.  
 
# PracTools 0.9

*  Errors in the help files for CVcalc2 and CVcalc3 corrected.  In CVcalc2, parameter 
descriptions for Wsq and Bsq were corrected. In CVcalc3, parameter description for delta1 was corrected. References to specific sections of Valliant, Dever, and Kreuter (2013) added for details on the formulas used in CVcalc2 and CVcalc3. Typo corrected in help file for NRadjClass. 
 
# PracTools 0.8

*   BW2stagePPSe changed to correct error in second term of between-PSU variance component estimate.  The formula used for the between component in v. 0.7 and earlier versions was: 

    \[\small {{v}_{PSU}}=\frac{1}{m\left( m-1 \right)}\sum\limits_{i\in s}{{{\left( \frac{{{{\widehat{t}}}_{i\pi }}}{{{p}_{i}}}-{{{\widehat{t}}}_{pwr}} \right)}^{2}}}-\sum\limits_{i\in s}{\frac{1-\pi _{i}^{*}}{{{\left( \pi _{i}^{*} \right)}^{2}}}{{{\widehat{V}}}_{i}}} \]
  
   where $\small \pi _{i}^{*}=m{{p}_{i}}$.  This estimator is somewhat biased when psu’s are selected with replacement.  An approximately unbiased estimator is  
  \[\small {{v}_{PSU}}=\frac{1}{m\left( m-1 \right)}\sum\limits_{i\in s}{{{\left( \frac{{{{\widehat{t}}}_{i\pi }}}{{{p}_{i}}}-{{{\widehat{t}}}_{pwr}} \right)}^{2}}}-\frac{1}{{{m}^{2}}}\sum\limits_{i\in s}{\frac{{{{\widehat{V}}}_{i}}}{p_{i}^{2}}}\]
   Since $\frac{1-\pi _{i}^{*}}{{{\left( \pi _{i}^{*} \right)}^{2}}}=\frac{1}{{{\left( m{{p}_{i}} \right)}^{2}}}-\frac{1}{m{{p}_{i}}}$, these two choices will be nearly equal when $\small {1}/{m{{p}_{i}}}\;$ is negligible compared to $\small {1}/{{{\left( m{{p}_{i}} \right)}^{2}}}\;$.  This will occur when all values of $m{{p}_{i}}$ are small.  The corrected estimator of $\small B^2$ is 
  \[\small {{\widehat{B}}^{2}}=\frac{m{{v}_{psu}}}{\widehat{t}_{pwr}^{2}} \]
 
# PracTools 0.7
*    wtdvar modified to handle NAs. 
    
*    mibrfss file added. File was used in  

    Valliant, R., and Dever, J. (2011), “Estimating Propensity Adjustments for Volunteer Web Surveys,” Sociological Methods and Research, 40, 105-137.  

# PracTools 0.6

*    Recompiled program with R CMD check practools per B. Ripley email of 2017/07/17. This corrects some disallowed CRLF endings in some lines.  Problem was old file called cleanup that was unnecessary.
    
# PracTools 0.5

  * BW2stagePPSe 
  
    Fixed help file for BW2stagePPSe: wrong description of pp parameter 

    Added k to output.  See Section 9.2.3 in Practical Tools book. 

  * BW2stagePPSe
  
    Added k1 and k2 to output.  See Section 9.2.4 in Practical Tools book. 

  **Various functions** 

*  Fixed if statements to eliminate warning when logic test is applied to combination of conditions. 
      
  * **Design effect functions**
  
*    All design effects can be requested in the deff function.  deff calls one of these depending on the value of a type parameter: deffCR, deffH, deffK, or deffS 
  * deffCR 
  
    Added Chen-Rust design effect for stratified, two-stage sampling with varying weights. See 

     Chen, S. and Rust, K. (2017). An Extension of Kish’s Formula for Design Effects to Two- and Three-Stage Designs with Stratification. \emph{Journal of Survey Statistics and Methodology}, 5(2), 111-130. 
  
# PracTools 0.4

*  clusOpt2. Calculation of CV when cal.sw = 1 corrected to be $\small \frac{{\tilde{V}}}{m\bar{n}}k\left[ 1+\delta \left( \bar{n}-1 \right) \right]$. $\small k$ was omitted in previous versions.
  
*  clusOpt3. Calculation of CV when cal.sw=1 or 2 corrected to be 
     \[\small \frac{{\tilde{V}}}{m\bar{n}\,\bar{\bar{q}}}\left\{ {{k}_{1}}{{\delta }_{1}}\bar{n}\,\bar{\bar{q}}+{{k}_{2}}\left[ 1+{{\delta }_{2}}\left( \bar{\bar{q}}-1 \right) \right] \right\}\].  
     ${{k}_{1}}$ and ${{k}_{2}}$ were omitted in previous versions.
     
* nhis, nhis.large, nhispart. Help files updated to correct descriptions of stratum variable. 
   
# PracTools 0.3

*    NRFUopt. Output for sample sizes in NRFUopt changed to unrounded numbers.  Note added to help file saying that unrounded sample sizes are used to compute total expected cost. 
BW2stagePPS, BW3stagePPS 

*    Description of pp in help file changed to “vector of one-draw probabilities for the PSUs; length is number of PSUs in population.”  The examples for these functions ran correctly in previous versions and use values of pp that conform to the corrected description in the help file. 

# PracTools 0.2

*  Classes for dub and strAlloc changed to power.htest.  

# PracTools 0.1

  * deffk
  
    New function added to compute the Kish design effect for unequal weighting.  
  
  * deffH 
  
    New function added to compute the Henry design effect that accounts for any gains in efficiency when estimating the total of a y variable using a GREG estimator based on a vector of covariates. 
  
  * deffS 
  
    New function added to compute the Spencer design effect that accounts for any gain in efficiency of estimating a total when a y variable depends on selection probabilities. 
  
  * nCont 
  
    Error trap added: alpha cannot be a vector. 
  
  * nContMoe 
  
    New function added to compute sample sizes in srswor for continuous variables based on margins of error.  
  
  * nhispart  
  
    Population added to package. An extract from NHIS 2003 used in an exercise in chapter 14 of Practical Tools for Designing and Weighting Survey Samples.  
  
  * nPropMoe 
  
    Help description corrected. Function only computes sample sizes for proportions. 
Error trap added to check that $\small 0<e<1$. 
  
  * NRadjClass 
  
    New function added to compute separate nonresponse adjustments in a set of classes.  Five methods for computing adjustments are included.  
  
  * pclass 
  
    New function added to fit a binary regression model for response probabilities and divide units into a specified number of classes. 
  
  * strAlloc 
  
    Warning now given if $\small {{n}_{h}}\ge {{N}_{h}}$. 
  
    Warning added: ch ignored if allocation = “neyman” 
  
  * Anticipated SE added for all allocations other than proportional. 
  
  * wtdvar 
  
    Warning now given if any weight is less than or equal to 0. 
  
    $\small {n}/{\left( n-1 \right)}$ factor added to agree with unbiased variance estimate in *srs*. 
     
# PracTools 0.0-2

  * BW2stagePPS: Description of pp parameter corrected in help file. 
  
  * BW3stagePPSe: Problem in example corrected.  Counts of SSUs per PSU in Ni corrected. 
  
  * strAlloc  
  
    Error message added if n.tot omitted when allocation is neyman. 
  
    If cost is non-NULL and alloc= "neyman", the cost parameter is ignored if used. Warning message added. 
  
    If alloc= "prop", n.tot is required.  Error trap added. 
  
    Check added to determine whether sum(ch) <= tot.cost.  Error message printed if condition violated. 
    
    n.tot should not be specified if the allocation is totvar or totcost.  An error message is issued if condition violated. 

  * nDomain 
  
    New function added to compute a simple random sample size using either a target coefficient of variation, $\small C{{V}_{0}}\left( d \right)$, or target variance, $\small {{V}_{0}}\left( d \right)$, for an estimated mean or total for a domain. 
  
  * CVcalc2 
  
    New function added to compute the coefficient of variation of an estimated total in a twostage design. Primary sampling units (PSUs) can be selected either with probability proportional to size (*pps*) or with equal probability. Elements are selected via simple random sampling (*srs*). 
  
  * CVcalc3 
  
    New function added to compute the coefficient of variation of an estimated total for a three-stage sample. PSUs can be selected either with varying probabilities and with replacement or with equal probabilities and with replacement. SSUs and elements within SSUs are selected by simple random sampling. The CV formula is appropriate for approximating the relvariance of the probability-with-replacement (*pwr*)-estimator of a total when the same number of SSUs is selected in each PSU and the same number of elements is selected within each sample SSU. 
     
