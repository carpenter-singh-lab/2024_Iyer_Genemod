# Manuscripts Data and Figures

Follow these steps to set up and run the analysis:

## Prerequisites

1. Clone and link the dependency repository:

```bash
# Clone the pilot data repository
cd ~/Desktop
git clone git@github.com:jump-cellpainting/pilot-cpjump1-data.git
cd pilot-cpjump1-data
git checkout v1.0.0

# Create symbolic link in analysis directory
cd /path/to/this/repo/manuscript_analysis/1.analysis
ln -s ~/Desktop/pilot-cpjump1-data .
```

## Environment Setup 

### R Environment
This project uses RStudio and renv for reproducibility. To get started:

1. Open the project by loading `genemod.Rproj` in RStudio
2. When prompted, allow renv to bootstrap (version 0.13.1)
3. Run `renv::restore()` to install required packages

### Python Environment for Supplementary Analysis
Before running the supplementary analysis, set up the Python environment:

```bash
# Navigate to supplementary directory
cd 5.supplementary

# Create environment from specification file
mamba env create -f environment.yml

# Activate the environment
mamba activate genemod_supplementary

# If environment.yml has been updated, run:
mamba env update -f environment.yml
```

## Running the Analysis

Execute the R notebooks in sequence:

```R
source("0.inspect-metadata/0.knit-notebooks.R", chdir = T)
source("2.inspect-analysis/0.knit-notebooks.R", chdir = T)
source("4.figures/0.knit-notebooks.R", chdir = T)
```

Then run the supplementary analysis:

```bash
cd 5.supplementary
bash run_analysis.sh
```