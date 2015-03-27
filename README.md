This repository contains the data and code necessary to reproduce the analysis in the manuscript "Spectating is stressful: observing a brief raptor attack elicits a glucocorticoid response in a social passerine" submitted to the Proceedings of the Royal Society B.  

Authors (affiliation):
- Blake C. Jones (Department of Biological Sciences, University of Memphis)
- Adam D. Smith (Department of Natural Resources Science, University of Rhode Island)
- Sara E. Bebus (Department of Biological Sciences, University of Memphis)
- Stephan J. Schoech (Department of Biological Sciences, University of Memphis)

To run the analysis:

1. Fork this repository
2. Open [R](http://www.r-project.org) on a computer with internet access
3. Install the [devtools](http://cran.r-project.org/package=devtools) and [knitr](http://cran.r-project.org/package=knitr) packages, if necessary 
  - install.packages("devtools")
  - install.packages("knitr")
4. Knit the file to produce an html document
  - knitr::knit2html('PATH/TO/Analysis_summary.Rmd')
5. Open the resulting html file (Analysis_summary.html) in your favorite browser

Note that we performed degrees of freedom adjustments in SAS.  We provide the SAS script (./SAS/final_analysis.sas) that produces the adjusted tests.  We also include the relevant output of the script as two `csv` files is the SAS directory.