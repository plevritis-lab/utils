#!/bin/bash

SEGMENTATION_METHOD="<TODO>"
NODAL_STATUS="<TODO>"

DATA_DIRECTORY="<TODO>"
CELESTA_DIRECTORY="<TODO>"
SCRIPT_DIRECTORY="<TODO>"

CONDITION="<TODO>"

CLINICAL_DATA="$DATA_DIRECTORY/clinical_annotations/$CONDITION/<TODO>.csv"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/<TODO>.csv"
THRESHOLDS_DIRECTORY="$CELESTA_DIRECTORY/thresholds/$CONDITION/$NODAL_STATUS"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION/$NODAL_STATUS"
QUANTIFICATIONS_DIRECTORY="$IMAGE_DIRECTORY/quantifications/$SEGMENTATION_METHOD"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"

THRESHOLD_GENERATOR_SCRIPT="$SCRIPT_DIRECTORY/generate_thresholds.py"
CELESTA_SCRIPT="$SCRIPT_DIRECTORY/apply_celesta.R"
VISUALIZATION_SCRIPT="$SCRIPT_DIRECTORY/visualize_dynamic_overlays.py"

python3 "$THRESHOLD_GENERATOR_SCRIPT" \
    --clinical_data "$CLINICAL_DATA" \
    --save_path "$THRESHOLDS_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX"

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
        SAMPLE_PROTEOMIC_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"
        SAMPLE_HISTOLOGY_DATA="$IMAGE_DIRECTORY/histology/${SAMPLE_NAME}.tif"
        SAMPLE_ASSIGNMENTS="$IMAGE_DIRECTORY/assignments/${SAMPLE_NAME}_assignments.csv"

        SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/<TODO>/image_1_seg.npy"

        open -a "Microsoft Excel" "$THRESHOLDS_DIRECTORY/${SAMPLE_NAME}_thresholds.csv"
        
        python3 "$VISUALIZATION_SCRIPT" \
            --assignments_path "$SAMPLE_ASSIGNMENTS" \
            --colormap_path "$SCRIPT_DIRECTORY/<TODO>.json" \
            --image_path "$SAMPLE_PROTEOMIC_DATA" \
            --histology_path "$SAMPLE_HISTOLOGY_DATA" \
            --mask_path "$SAMPLE_SEGMENTATION" \
            --use_$SEGMENTATION_METHOD
    fi
done