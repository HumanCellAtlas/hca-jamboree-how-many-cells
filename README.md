# Prospective experimental design task : How Many Cells? 

The purpose of this task is to help answer the question that we are always asked (or ask ourselves!): How many cells do I need to sequence for my experiment?

We know that the ability to separate cells into clusters is dependent on many factors, including:

* Total cell number in the experiment
* Frequency of the cluster
* Sequencing depth
* Transcriptional 'distinctness' of the clusters.

Your task is to build a model to predict the separability of any two pairs of clusters, for any dataset. You can see exact specifications for inputs and outputs in the task README. This repository contains example scripts for the HCA jamboree problem on cluster separability. Referenced data can be found on the Jamboree server under `/data/tasks/how_many_cells`

#### Subtask 1: An estimator of cluster separability

 * Think of different definitions of cluster 'separability'. How can we design a metric that has minimal assumptions on the data, and is reflective of the power in the experimental design? We provide a reference implementation (essentially, the fraction of KNN in the 'correct' cluster), but welcome additional ideas and suggestions.
 
#### Subtask 2: Generate and analyze separability, as a function of cell number and sequencing depth

 * Generate separability curves for pairs of clusters, across a variety of cell depths. We provide a vignette to demonstrate this on a human pancreas dataset, but encourage you to perform this task across multiple datasets.
 * Intuitively, what are key determinants of separability across datasets, depths, and pairs of clusters?

#### Subtask 3: Predict separability curves, for a new experiment

 * Based on a subset of cells drawn from two clusters, and their frequencies, predict the separability as a function of cell count.
 * Consider linear and non-linear predictors. Which features can be extracted from the sample cell population to make this informative? 
 * Use our definition of separability and beat the reference implementation. 

#### Subtask 4: Identify a minimal, and interpretable, feature set to guide future experimental design.

 * The model in T3 can use a wide range of features that can be derived from the cell populations. An important question is to condense this information into a minimal set of interpretable features. 
 * Try to define 2-3 quantities that really matter. To be considered (for example): frequency, log fold change difference of DE genes, #DE genes. 
 * The ultimate aim is to define effect size and frequency-like measures similar to traditional power calculation designs. 

### Available, pre-clustered, datasets
 * ~26,000 mouse retinal bipolar cells (Shekhar et al., Cell, 2016)
 * ~10,000 human pancreatic islet cells (Baron et al, Cell Systems 2016)
 * ~33,000 human PBMCs (Zheng et al, Nature Communications, 2017)
 * Larger datasets to appear on Day 2

### First Deliverable : Calculate actual separability for two clusters in a dataset

**Input:** Gene expression matrix for a full dataset, with pre-determined clustering

**Output:** Separability of any two clusters

**Input arguments to script:**

* Loom file of entire dataset, contains cluster ID for each cell
* Names of the two clusters to determine separability (as separate values)

**Output:**

* A single number: the "separability" of the two clusters, output to console

Example: `Rscript --vanilla actual_separability.R human_pancreas.loom activated_stellate quiescent_stellate`


### Second Deliverable : Predict separability of two clusters, as a function of total cell number in the dataset

**Input:** Gene expression only for example cells that are members of the two clusters, and their cluster frequencies in the full dataset

**Output:** Predicted separability of the clusters as a function of total cell number.

**Input arguments to script:**

* Loom file containing example cells representing the two clusters. Contains the cluster ID for each cell
* Cluster IDs and frequencies for the two clusters

**Output:** Predicted separability of the clusters as a function of total cell number.
  
* Frequency of the two clusters, in the total dataset (as separate values)
* Path to a TSV file, which will contain two columns.
  1. Number of cells (ranging from 1,000 to 100,000 - with a step of 1,000)
  2. Predicted separability

Example: `Rscript --vanilla predicted_separability.R human_pancreas_activated_quiescent_stellate.loom activated_stellate 0.03 quiescent_stellate 0.02 predictions.tsv`

### Additional notes:

* All data files are stored according to the loom specification, an HDF5-based format. Each contains a UMI matrix, cell_names and cluster ID as column attributes, and gene_names as a row attribute.
* For our example "predicted separability" model, we construct a dummy function that considers only the number of DE genes (>2-fold change) as input.

### Example task execution

* See [hcajamboree_howmanycells_vignette.html](https://cdn.rawgit.com/HumanCellAtlas/hca-jamboree-how-many-cells/ajc-upload/hcajamboree_howmanycells_vignette.html) for a worked example of this task. 
