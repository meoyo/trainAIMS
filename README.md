Fist of all, thank you for your interest in AIMS.

I haven't created a R package per se for the R functions necessary to generate the AIMS' like models, but I have a couple of R files that could be used to generate such models. I wrapped them up and I included them in this folder.

Before you start, you need to install those R packages :

1) e1071  
2) gplots  
3) ROCR  
4) multicore  

that are necessary to run the scripts. It involves mostly running install.packages(package) from within R.

You would also need the Rgtsp package from Vlad Popovici :

1) Use the Rgtsp_0.0-5.tar.gz included in the package
2) launch in the terminal : R CMD INSTALL Rgtsp_0.0-5.tar.gz

There is a lot of stuff and function that are not necessary for you in the folder.

In principle, you will only need to source("trainAIMS.R") and run the trainAIMS function. Please just have a quick look at the function in the file trainAIMS.R. Hopefully you will be able to understand how to launch the function by reading the comments.

Just let me know if there is something you don't understand (eric.r.paquet@gmail.com).

The way we are generating our AIMs' models is somehow unique (at least in comparison to other top scoring pairs models) and we haven't published this stuff yet. Consequently, I will appreciate if you could cite :

Absolute assignment of breast cancer intrinsic molecular subtype.  
Paquet ER, Hallett MT.  
J Natl Cancer Inst. 2014 Dec 4;107(1):357. doi: 10.1093/jnci/dju357. Print 2015 Jan.  
PMID: 25479802  

in case you publish models generated using this package.

It is a good idea to call the function set.seed(your_favorite_number_here) before launching the trainAIMS function so your cross-validation splits are consistent from one run to another.

Whatever results you are getting with your project, I will be interested to know what type of results you obtained. Thanks!

Have a good day  
Sincerely,  
Eric (eric.r.paquet@gmail.com)  
