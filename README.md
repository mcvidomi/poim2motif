# poim2motif
Assessing motifs in Positional Oligomer Importance Matrices (POIMs).

# required installations

1. shogun toolbox<br />
    http://www.shogun-toolbox.org/doc/en/3.0.0/installation.html<br />
    install with python interface
2. openopt<br />
    http://openopt.org/Install
3. R
    https://cran.r-project.org/mirrors.html <br />
    install the additional package bioconductor package by enter in R<br />
        source("http://bioconductor.org/biocLite.R")<br />
        biocLite("seqLogo")<br />


# tutorial

"run_toy.py" is the toy Example
"run_real1.py" is the real Example. For this real experiment you need to download one folder of the real data from: http://www.fml.tuebingen.mpg.de/raetsch/projects/lsmkl and adjust the datapath in "run_real1.py"


1. create POIM<br />
    by set CURRENT_TASK to TASK_1

    
2. compute Motif<br />
    by set CURRENT_TASK to TASK_2

With "path" you specify the folder path, including all results.
You can choose the Name of the Experiment by experiment_name, which will be the name of the folder including all results.

If you have any questions, please write me a mail at marina.vidovic@tu-berlin.de.

Good luck!
