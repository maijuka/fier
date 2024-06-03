# Fier instructions

This is the code for paper "Fast Redescription Mining Using Locality-Sensitive Hashing". Repository includes code from the ReReMi package available at https://pypi.org/project/python-clired/.

## Setting up the environment 

Create a pipenv using the given Pipfile, or install Python version 3.6 or higher. You should also install at least scipy, numpy, matplotlib, scikit-learn, wxPython and cartopy.

## Source code

The algorithms presented in the paper are implemented in the file classMinerLSH.py. 

The file used to run the algorithm is called exec_clired.py.

You do not need to edit any of the code files, but to run the algorithm you need to use a settings file to define the data and settings used.

## How to set the input data

The same settings files can be used for any dataset. The dataset used is set with the following parameters 

* LHS_data, insert the file name of the left hand side data
* RHS_data, same as above but for the right hand side data

The names of the datasets match the ones in the paper, but make sure you have the correct paths for them.

To save the results to a folder, add a path to the parameter

* result_rep

Some experiments also use pre-computed initial pairs to extend. If there is a parameter called 

* pairs_store

Make sure that you have the correct path for the file given there.

Also make sure that the path to the file in 

* fields_rdefs

is correct.

## How to set the parameters to run Fier

Following parameters are for all datasets except DentalW and DentalA (their parameters are specified in the paper).

### The parameters that were the same for all experiment runs (also ReReMi and ReReMiBKT) are

* max_var_s0 = 0
* max_var_s1 = 1
* min_itm_out = 0.3
* min_itm_in = 0.1
* max_inits = 100
* init_minscore = 0

### To use LSH this parameter needs to be set to 1

* lsh_extensions = 1

### With LSH the following parameters are the same for all experiments

* lsh_min_itm_out = 1.5
* method_buckets = similar-height
* nb_ext_buckets = 2
* lsh_bucket_mult = 5
* lsh_nb_total_exts = 10000

### For Fier_init parameters b_j and r_j set (sections 3.2 and 3.4)

* lsh_nb_bands = 40
* lsh_nb_rows = 10

### For Fier_init with numerical data use also

* method_buckets = similar-height
* nb_buckets = 40

### For Fier_ext parameters b_h and r_h set (section 3.3 and 3.4)

* lsh_nb_bands_ext = b_h
* lsh_nb_rows_ext = r_h

### The experiments comparing LSH extension quality with ReReMi (section 3.3, second paragraph) have also

* ext_once = 1

Which means each initial pair will be extended only once at most.

## How to set the parameters to run ReReMi and ReReMiBKT

If you want to run ReReMi, make sure you don't set 
### To run ReReMiBKT for the initial pairs, set 

* method_buckets = similar-height
* nb_buckets = 40

### If you want to only mine pairs with ReReMi or ReReMiBKT set

* only_pairs = 1

## Running the experiments

With a settings file "settings.xml", you can run any of the experiments with

python exec_clired.py settings.xml

The demo.zip folder includes test data (LHS and RHS for left hand side and right hand side data tables) with pre-mined pairs ( _pairs.txt). The files are named according to the datasets used in the paper. Files with no dataset name (fields_rdefs_custom.xml) are used for all experiments. 
