#!/bin/zsh

for dir in ../assembloids/*/; do
    sample_name=$(basename "$dir")

    python3 segmentation/quantify_expression.py \
        --image_path ../assembloids/$sample_name/data/$sample_name.tif \
        --mask_path ../assembloids/$sample_name/full/mesmer/"Vim + EpCAM"/image_1_seg.npy \
        --panel_path ../assembloids/channel_names.txt --use_mesmer
done