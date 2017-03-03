# InitExp : Play and visualize EM Algorithm on Gaussian Mixture datasets

This R package is a collection of wrapper functions and test scripts that investigate the property of the Expectation Maximization algorithm. Specifically, it demonstrates some intriguing behavior of EM over a highly non-convex likelihood function.

#Install
To install directly from github, open a terminal, type R, then

    devtools::install_github('htso/InitExp')

#Dependencies
You need my Hext package from github, from a terminal, type 

    devtools::install_github('htso/Hext')

as well as these packages on CRAN,    

    install.packages("Hext", "mvtnorm", "corpcor", "ellipse")

#Datasets
I include a couple of generic datasets, which are different flavor of gaussian mixture from low to high dimension. To load a dataset, just type

    data(simdat2single)

#Run
I provide twp demos in the /demo subfolder. To run a demo, 

    demo("init-demo", package="InitExp")

#WARNING
Some of these scripts may take a long time to finish.

#Platforms
Tested it on Linux (ubuntu 14.04).





