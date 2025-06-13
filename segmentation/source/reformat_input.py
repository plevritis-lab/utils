import argparse
import os
import shutil
from tifffile import TiffFile, imread, imwrite

def condense_channels(image_path, panel_path, channels_to_remove):
    """condenses the channels of a multiplexed image by removing specified ones while updating the corresponding panel file

    args:
        image_path (str): path to image file in .tif / .tiff format; \
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

def reformat_images(image_path):
    """reformats an image by transposing its dimensions, saving to a new directory
    
    args:
        image_path (str): path to image file in .tif / .tiff format; \
                          image should be of shape (y, x, c) and will be reformatted to (c, y, x)
    """
    
    image = imread(image_path)

    image = image.transpose(2, 0, 1) # (c, y, x)

    sample_name = os.path.basename(image_path)
    save_directory = os.path.join(os.path.dirname(image_path), "reformatted_images")
    save_path = os.path.join(save_directory, sample_name)

    os.makedirs(save_directory, exist_ok = True)
    imwrite(save_path, image, check_contrast = False)

def prepare_input(data_directory):
    """prepares the input directory by organizing image files into subdirectories

    args:
        data_directory (str): the path to the directory containing the image files
    """

    for image_name in os.listdir(data_directory):
        if image_name.endswith((".tiff", ".tif", ".qptiff")):
            image_path = os.path.join(data_directory, image_name)
            sample_name = os.path.splitext(image_name)[0]
            sample_directory = os.path.join(data_directory, sample_name, "data")

            os.makedirs(sample_directory, exist_ok = True)

            shutil.move(image_path, os.path.join(sample_directory, image_name))




def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for formatting input directories prior to applying mesmer or cellpose")
        
    parser.add_argument("-d", "--data_directory", type = str, help = "path to the directory containing the images to be segmented; \
                                                                      images should be in the .tiff, .tif, or .qptiff format")  
    
    return parser.parse_args()


import tifffile
import numpy as np

def reformat_image(image_path):
    """reformats an image by transposing its dimensions, overwriting the original file
    
    args:
        image_path (str): path to image file in .tif / .tiff format; \
                          image axes must contain x, y, c, and optionally z in the metadata; \
                          image will be reformatted to (c, y, x) or (z * c, y, x) if z is present
    """

    with TiffFile(image_path) as tif:
        image = tif.asarray()
        axes = tif.series[0].axes

    print(axes)
    print(image.shape)

    assert(all(a in axes for a in ("C", "X", "Y")))

    if "Z" in axes:
        image = np.moveaxis(image, [axes.index("Z"),
                                    axes.index("C"),
                                    axes.index("Y"),
                                    axes.index("X")],
                            [0, 1, 2, 3])
        image = image.reshape(-1, image.shape[2], image.shape[3])
    
    else:
        image = np.moveaxis(image, [axes.index("C"),
                                    axes.index("Y"),
                                    axes.index("X")],
                            [0, 1, 2])
    
    imwrite(image_path, image)
    
def main():
    arguments = parse_arguments()

    # prepare_input(arguments.data_directory)

    reformat_image("/Users/rohit/Downloads/example/521_S1_reg024.tif")

main()