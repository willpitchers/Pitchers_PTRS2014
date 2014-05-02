Pitchers_PTRS2014
=================

This is the repository for the `R` Code and Data for Pitchers et. al. 2014 soon to be published (July 2014) a special issue (Evolutionary Quantitative Genetics) of Philosophical Transactions of the Royal Society. We currently have a pre-print of a version of the manuscript [here on biorxiv](http://dx.doi.org/10.1101/002683). In addition to all of the code and data being here allowing you to reproduce all of the results from the paper (up to the limits of MCMC anyways), it also contains the database we generated on genetic covariance and correlation matrices that we hope will be of broad use to the Evolutionary Quantitative Genetics community. 

Here is a short guide to navigate this folder:

1. For the analyses of evolutionary potential using measures extracted from empirical **G** matrices, you can [review](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/gmax_analysis.md) or reproduce our analyses. You can also download the [gmax_analysis.R](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/gmax_analysis.R) file to rerun it on your local machine. This folder also contains the [gmax_analysis.Rmd](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/gmax_analysis.Rmd) file in case you prefer starting with the markdown + `R` code and want to use [knitr](https://github.com/yihui/knitr).

2. If you are interested in the analyses we performed examining rates of evolution or estimates of selection, thse can also be found in the script folder, with the overview of the [analysis and results here in the markdown file](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/selection_%26_rate_analyses_final.md). The [`.R`](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/selection_%26_rate_analyses_final.R) file with just the code, and the `.Rmd` containing both the `.R` and  `.md` markdown can be found [here](https://github.com/DworkinLab/Pitchers_PTRS2014/blob/master/Scripts/selection_%26_rate_analyses_final.Rmd).

3. If you just interested in the data for your own analyses, please find it in the [data subfolder](https://github.com/willpitchers/Pitchers_PTRS2014/tree/master/Data). Inside you will see additional sub-folders which contain small text files for each matrix as spelled out in the name.md file mentioned above. You can also use some of the `R` code we have set-up above, to extract the databases directly into R with the `.R` objects in that folder.

