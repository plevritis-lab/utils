#!/bin/bash

SEGMENTATION_METHOD="mesmer"

CONDITION="normal"

DATA_DIRECTORY="$HOME/Desktop/attardi_analysis/data"
CELESTA_DIRECTORY="$HOME/Desktop/attardi_analysis/celesta"
SCRIPT_DIRECTORY="$HOME/Desktop/utils/celesta/source"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/normal_signature_matrix.csv"
THRESHOLDS_DIRECTORY="$CELESTA_DIRECTORY/thresholds/$CONDITION"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION"
QUANTIFICATIONS_DIRECTORY="$IMAGE_DIRECTORY/quantifications/$SEGMENTATION_METHOD"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"

THRESHOLD_GENERATOR_SCRIPT="$SCRIPT_DIRECTORY/generate_thresholds.py"
CELESTA_SCRIPT="$SCRIPT_DIRECTORY/apply_celesta.R"
VISUALIZATION_SCRIPT="$SCRIPT_DIRECTORY/visualize_dynamic_overlays.py"

python3 "$THRESHOLD_GENERATOR_SCRIPT" \
    --image_directory "$IMAGE_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX" \
    --save_path "$THRESHOLDS_DIRECTORY"

Rscript "$CELESTA_SCRIPT" \
    --data_directory "$QUANTIFICATIONS_DIRECTORY" \
    --save_path "$ASSIGNMENTS_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX" \
    --thresholds_directory "$THRESHOLDS_DIRECTORY" \
        > /dev/null 2> "$SCRIPT_DIRECTORY/celesta_warnings.txt"

open -a "Microsoft Excel" "$SIGNATURE_MATRIX"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        SAMPLE_PROTEOMIC_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.qptiff"
        SAMPLE_ASSIGNMENTS="$IMAGE_DIRECTORY/assignments/${SAMPLE_NAME}_assignments.csv"

        SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/E_CADHERIN/image_1_seg.npy"

        open -a "Microsoft Excel" "$THRESHOLDS_DIRECTORY/${SAMPLE_NAME}_thresholds.csv"
        
        python3 "$VISUALIZATION_SCRIPT" \
            --assignments_path "$SAMPLE_ASSIGNMENTS" \
            --colormap_path "$SCRIPT_DIRECTORY/colormaps/colormap_attardi.json" \
            --image_path "$SAMPLE_PROTEOMIC_DATA"
    fi
done