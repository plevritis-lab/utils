# segmentation

the segmentation module will process an entire directory of TIF / TIFF / QPTIFF files using either mesmer or cellpose; please note that segmentation is expensive - therefore, it is recommended to run the following pipeline on an HPC cluster like sherlock

*prelude*: note that cellpose, but not mesmer, works natively on apple's M chips

source files:

1. `source/format_directories.py`: organizes a directory of TIF / TIFF / QPTIFF images into nested subdirectories
2. `source/apply_segmentation.py`: runs either cellpose or mesmer on proteomic images; typically, cellpose should use cytoplasmic + nuclear markers and mesmer should use membrane + nuclear markers; the user has the flexibility to specify multiple such markers when selecting either algorithm
3. `source/quantify_expressions.py`: using the segmentation masks generated, quantifies whole-cell measurements, including morphological structure and fluorescent intensity

bash files:

1. `bash/segmentation_template.sh`: batch processes segmentation across multiple samples 
2. `bash/quantify_expressions_template.sh`: batch processes measurement quantification across multiple samples