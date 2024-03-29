# NoLabelPFA

Here are the preparation steps for the no label pfa

# Preparation in R (PrepareSeurat2PFAinput.R)

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

# No Label PFA Preparation:

copy the scripts 
* **01_File_Preparation_Script.ipynb** (in No_label_PFA/01_File_Preparation_Script.ipynb)
* **02_No_label_PFA.ipynb** (in No_label_PFA/02_No_label_PFA.ipynb)

and the unzipped **src** directory (available as zip file in No_label_PFA)  
as well as the previously generated PFA_input table  
into the same directory 

# Preparation in Python (No_label_PFA/01_File_Preparation_Script.ipynb) 

This step allows creating a test file, which is an optional step that is not necessary for the full analysis but can be performed as a quick check.
Alternatively, instead of creating a new testfile, you can download and unzip our testfile which was generated using the result of the R preparation of the data published by Sole-Boldo et al. (2020) (https://www.nature.com/articles/s42003-020-0922-4)

Since the file name of the R preparation depends on your input data and your variable names, you have to  
-> enter your file name in `df = pd.read_csv( )` 

This script performs a full analysis.  
**If you would like to perform a full analysis**, you have to run all steps of the **01_File_Preparation_Script.ipynb** script, *except* for the second step of the script (#creating TESTFILE: FOR TESTING PURPOSES ONLY!!! PLEASE DO NOT USE IF YOU WANT TO PERFORM A FULL ANALYSIS) 

**If you would like to create a test file instead**, you will have to perform the second step of the script (#creating TESTFILE: FOR TESTING PURPOSES ONLY!!! PLEASE DO NOT USE IF YOU WANT TO PERFORM A FULL ANALYSIS) and change the code as described below.  



## OPTIONAL use the provided testfile
copy the unzipped file (Test_File.csv) into the same directory as 01_File_Preparation_Script.ipynb
and change the code accordingly:

### in #creating preprocessed_data.csv
use:
`filename = 'Test_File.csv'` 

instead of:
`filename = 'Young_as_Young_and_Old_as_Old.csv'` 

--> this will create **preprocessed_data.csv**

### in #creating comparison_labels.csv
use:
`data = pd.read_csv('Test_File.csv', nrows=1, header=None)`

instead of:
`data = pd.read_csv('Young_as_Young_and_Old_as_Old.csv', nrows=1, header=None)`

--> this will create **comparison_labels.csv**

### in #creating Genname.csv
use:
`input_file = "Test_File.csv"`

instead of:
`input_file = "Young_as_Young_and_Old_as_Old.csv"`

--> this will create **Genname.csv**

# No Label PFA (No_label_PFA/02_No_label_PFA.ipynb) 

For the second step, the directory **src** (available as zip file) needs to be in the same directory as the input files that have been created with 01_File_Preparation_Script.ipynb.

Since **preprocessed_data.csv** was either generated for the full analysis or for the test file (if you performed step 2 of 01_File_Preparation_Script.ipynb and changed the code accordingly) the subsequent steps require no further adjustment of the input name.  
Since you can choose the clusters that best suit your analysis, you need to adujst the code in the steps **Analyzing Cluster Differences** and **Get Mutual Information**, and the **Validation** as well as in the **Shaply Explanation** and/or the **Tree Explanation** steps

For instance, for cluster 0 and 2: **find_cluster_differences(path_original_data=path_original_data,clusters=[0,2])** and **get_mutual_information(path_original_data,clusters=[0,2])**, and the **validate_feature_selection(path_original_data,number_sweeps=20, n_highest_mutual_information=5, feature_selection=0, clusters=[0,2])** as well as in the **shaply_explanation(path_original_data, n_highest_mutual_information=1, clusters=[0,2])** and/or the **tree_explanation(path_original_data, n_highest_mutual_information=1, clusters=[0,2], min_samples_leaf=100)** steps

## OPTIONAL
You can choose **tsne** instead of **umap**  
and **hdbscan** instead of **dbscan**  
by changing the code of **02_No_label_PFA.ipynb** accordingly

## Analyzing Cluster Differences
in `compare_dbscan_labels("comparison_labels.csv")` (or `compare_hdbscan_labels("comparison_labels.csv")` if you opted for HDBSCAN)  
you will receive a bar plot with each bar visualizing the cells of a cluster.  
-> chose 2 clusters you would like to compare (for instance cluster 0 and cluster 2):  
`find_cluster_differences(path_original_data=path_original_data,clusters=[0,2])`

If you receive an error warning, try selecting bigger clusters (taller bars in the bar plot)

## Get Mutual Information
chance the clusters accoring to your selected clusters (e.g., cluster 0 and cluster 2):  
`get_mutual_information(path_original_data,clusters=[0,2])`

## Validation (# feature_selection: PFA = 0, random features = 1 or all features = 2)
change the clusters according to your selected clusters (e.g., cluster 0 and 2)
`validate_feature_selection(path_original_data,number_sweeps=20, n_highest_mutual_information=5, feature_selection=0, clusters=[0,2])`

## Shaply Explanation
change the clusters according to your selected clusters (e.g., cluster 0 and 2)
`shaply_explanation(path_original_data, n_highest_mutual_information=1, clusters=[0,2])`

## Tree Explanation
change the clusters according to your selected clusters (e.g., cluster 0 and 2)
`tree_explanation(path_original_data, n_highest_mutual_information=1, clusters=[0,2], min_samples_leaf=100)`
