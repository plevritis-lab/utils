#!/bin/bash

SEGMENTATION_METHOD="mesmer"

CONDITION="primary_tumors/edges/node_positive"

DATA_DIRECTORY="$HOME/Desktop/plevritis_analysis/data"
CELESTA_DIRECTORY="$HOME/Desktop/plevritis_analysis/celesta"
SCRIPT_DIRECTORY="$HOME/Desktop/EPOCHS/celesta/source"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/edge_positive_signature_matrix.csv"
THRESHOLDS_DIRECTORY="$CELESTA_DIRECTORY/thresholds/$CONDITION"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION"
QUANTIFICATIONS_DIRECTORY="$IMAGE_DIRECTORY/quantifications/$SEGMENTATION_METHOD"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"
VISUALIZATION_DIRECTORY="$IMAGE_DIRECTORY/visualizations"

THRESHOLD_GENERATOR_SCRIPT="$SCRIPT_DIRECTORY/generate_thresholds.py"
CELESTA_SCRIPT="$SCRIPT_DIRECTORY/apply_celesta.R"
DYNAMIC_VISUALIZATION_SCRIPT="$SCRIPT_DIRECTORY/visualize_dynamic_overlays.py"
STATIC_VISUALIZATION_SCRIPT="$SCRIPT_DIRECTORY/visualize_assignments.py"

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

# open -a "Microsoft Excel" "$SIGNATURE_MATRIX"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" && "$SAMPLE_NAME" != "visualizations" && "$SAMPLE_NAME" != "topology" ]]; then
        SAMPLE_PROTEOMIC_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"
        SAMPLE_HISTOLOGY_DATA="$IMAGE_DIRECTORY/histology/${SAMPLE_NAME}.tif"
        SAMPLE_ASSIGNMENTS="$IMAGE_DIRECTORY/assignments/${SAMPLE_NAME}_assignments.csv"

        SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/PDGFRB + CD3 + CYTOKERATIN/image_1_seg.npy"

        # open -a "Microsoft Excel" "$THRESHOLDS_DIRECTORY/${SAMPLE_NAME}_thresholds.csv"

        # python3 "$DYNAMIC_VISUALIZATION_SCRIPT" \
        #     --assignments_path "$SAMPLE_ASSIGNMENTS" \
        #     --colormap_path "$SCRIPT_DIRECTORY/colormaps/plevritis/colormap_primary_plevritis.json" \
        #     --image_path "$SAMPLE_PROTEOMIC_DATA" \
        #     --histology_path "$SAMPLE_HISTOLOGY_DATA" \
        #     --mask_path "$SAMPLE_SEGMENTATION" \
        #     --apply_$SEGMENTATION_METHOD

        python3 "$STATIC_VISUALIZATION_SCRIPT" \
            --assignments_path "$SAMPLE_ASSIGNMENTS" \
            --colormap_path "$SCRIPT_DIRECTORY/colormaps/plevritis/colormap_primary_plevritis.json" \
            --save_path "$VISUALIZATION_DIRECTORY"
    fi
done