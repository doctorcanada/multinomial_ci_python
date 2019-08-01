# multinomial_ci.py

This is a Python implementation of the Sison-Glaz (1995) method for 
constructing simultaneous confidence intervals for multinomial proportions. 
This is a translation of the source code used to implement the "MultinomialCI" 
R package (source code: https://github.com/cran/MultinomialCI, originally
prepared by Dr. Pablo Villacorta), which itself was an R-based implementation 
of a SAS package coded by May and Johnson (2000).
An example usage is provided in the definition of the "main" function at
the bottom of this script.

## References
Sison, C. &  Glaz, J. (1995). Simultaneous Confidence Intervals and 
   Sample Size Determination for Multinomial Proportions, Journal of the 
   American Statistical Association, 90:429, 366-369, 
   DOI: 10.1080/01621459.1995.10476521
   
May, W., & Johnson, W. (2000). Constructing two-sided simultaneous confidence 
   intervals for multinomial proportions for small counts in a large number of 
   cells. Journal of Statistical Software, 5:6, 1-24. 
   doi:http://dx.doi.org/10.18637/jss.v005.i06
   [Online: , see https://www.jstatsoft.org/article/view/v005i06]
   
R package for MultinomialCI, by P. Villacorta:
https://cran.r-project.org/package=MultinomialCI
   
PLEASE NOTE THAT THE R PACKAGE HAS A MINOR BUG which has been reported
to the original package author.

## Created by
Brian Canada, PhD (bcanada@uscb.edu) 

Date of original commit: 31 July 2019
