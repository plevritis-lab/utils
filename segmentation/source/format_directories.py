import argparse
import os
import shutil

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
    
def main():
    arguments = parse_arguments()

    prepare_input(arguments.data_directory)

main()