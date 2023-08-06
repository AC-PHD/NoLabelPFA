# NoLabelPFA

Here are the preparation steps for the no label pfa

# Preparation in R

Our example uses the data published by Sole-Boldo et al. (2020) (https://www.nature.com/articles/s42003-020-0922-4), 
the RDS-file is available via the GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130973

The script will create subsets according to condition
during the preparation, it will remove ribosomal and mitochondrial genes 
and create PFA input tables which are labeled with A) the condition name and B) the a number for the condition (see the chapter NAME variables)
if you do not intend to use numbers as labels you can ignore the chapter "Create Table using numbers as labels"

in the section "# NAME variables", you need to:
* set the variable names according to the conditions of your samples (in our example "young" and "old")
* set the variable label numbers (e.g., "0" and "1")
* read in the RDS file
* and set the criterion variables according to the conditions of the data set 
  (the column names in the metadata of the Seurat object e.g., "age" and "celltype.age"

-> this will result in:
* a PFA_input table using the conditions as label (recommended) 
* a PFA_input table using the label numbers as label

# Preparation in Python

