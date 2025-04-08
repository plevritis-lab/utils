#!/bin/bash

SEGMENTATION_METHOD="<TODO>"

PANEL_PATH="<TODO>"
IMAGE_DIRECTORY="<TODO>"
SEGMENTATION_SCRIPT="<TODO>"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.<TODO>"
    
    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$SEGMENTATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --panel_path $PANEL_PATH \
            --nuclear_channel "<TODO>" \
            --segment_channel "<TODO>" \
            --apply_$SEGMENTATION_METHOD
    fi
done