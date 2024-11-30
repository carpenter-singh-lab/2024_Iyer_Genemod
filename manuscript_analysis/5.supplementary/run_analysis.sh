#!/bin/bash

# # Create and activate conda environment
# echo "Setting up conda environment..."
# mamba env create -f environment.yml
# source activate genemod_supplementary

# Create output directory if it doesn't exist
mkdir -p output

echo "1. Processing compound predictions and merging with repurposing annotations..."
python 1.well_predicted_compounds.py

echo "2. Performing over-representation analysis (ORA) for compound-gene and gene-pathway relationships..."
python 2.genemod_ora.py

echo "3. Analyzing chemical similarity and diversity of compounds targeting specific genes..."
python 3.chemsim.py

echo "4. Generating enrichment analysis tables for genes and pathways..."
python 4.plot_ora.py

echo "5. Processing transformer predictions for compound formulations..."
python 5.get_full_predictions.py

echo "6. Combining configuration files from grid search experiments..."
python 6.concat_configs.py

echo "7. Merging learning logs from grid search experiments..."
python 7.concat_learning_log.py

echo "8. Analyzing and visualizing learning metrics across different hyperparameters..."
python 8.learning_summary.py

echo "Analysis complete! Results can be found in the output directory."