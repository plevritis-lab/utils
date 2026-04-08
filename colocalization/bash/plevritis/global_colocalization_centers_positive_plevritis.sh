#!/bin/bash

CONDITION="primary_tumors/centers/node_positive"

DATA_DIRECTORY="$HOME/Desktop/plevritis_analysis/data"
CELESTA_DIRECTORY="$HOME/Desktop/plevritis_analysis/celesta"
SCRIPT_DIRECTORY="$HOME/Desktop/EPOCHS/colocalization/source/global_colocalization"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/center_positive_signature_matrix.csv"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"
VISUALIZATION_DIRECTORY="$IMAGE_DIRECTORY/visualizations"
COLOCALIZATION_DIRECTORY="$IMAGE_DIRECTORY/colocalizations"

GLOBAL_COLOCALIZATION_SCRIPT="$SCRIPT_DIRECTORY/generate_condition_summaries.R"

Rscript "$GLOBAL_COLOCALIZATION_SCRIPT" \
    --data_directory "$ASSIGNMENTS_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX" \
    --save_path "$COLOCALIZATION_DIRECTORY"