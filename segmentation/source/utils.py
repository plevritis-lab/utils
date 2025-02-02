import numpy as np
import os
from skimage.io import imread, imsave
from tifffile import TiffFile
from xml.etree import ElementTree

def extract_proteomic_panel(image_path, panel_path = None):
    """extracts the protein panel used in CODEX experiments by parsing through the underlying QPTIFF metadata

    args:
        image_path (str): file path that points to the underlying location of the QPTIFF
        panel_path (str, optional): file path that points to the underlying location of the channel_names.txt file; \
                                    this argument should only be supplied when metadata parsing fails
    """

    protein_panel = []

    if panel_path:
        with open(panel_path, "r") as panel:
            protein_panel = [m.rstrip() for m in panel]

    else:
        with TiffFile(image_path) as qptiff:
            for page in qptiff.series[0].pages:
                image_metadata = page.tags["ImageDescription"].value
                image_metadata = ElementTree.fromstring(image_metadata)

                protein = image_metadata.find("Biomarker").text

                protein_panel.append(protein)

    return protein_panel

def condense_channels(image_path, panel_path, channels_to_remove):
    """condenses the channels of a multiplexed image by removing specified ones while updating the corresponding panel file

    args:
        image_path (str): path to image file in .tif / .tiff format; 
                          image should be of shape (c, y, x) or (z, y, x, c)
        panel_path (str): path to the panel file containing channel names
        channels_to_remove (list): list of channel names to be removed from the image and panel file
    """

    with open(panel_path, "r") as panel:
        protein_panel = [m.rstrip() for m in panel]

    image = imread(image_path)

    if image.shape[-1] == 4: # (z, y, x, c)
        image = image.transpose(0, 3, 1, 2).reshape(-1, image.shape[1], image.shape[2])

    for index in range(len(protein_panel) - 1, -1, -1):
        if protein_panel[index] in channels_to_remove:
            image = np.delete(image, index, axis = 0)
            del protein_panel[index]

    sample_name = os.path.basename(image_path)
    save_directory = os.path.join(os.path.dirname(image_path), 
                                    "condensed_images")
    save_path = os.path.join(save_directory, sample_name)
    
    os.makedirs(save_directory, exist_ok = True)
    imsave(save_path, image, check_contrast = False)

    condensed_panel_path = os.path.join(os.path.dirname(panel_path), "condensed_channel_names.txt")
    with open(condensed_panel_path, "w") as panel:
        panel.write("\n".join(protein_panel))