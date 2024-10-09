import cv2
import numpy as np
import os
import tifffile as tiff

def darken_histology(image_path, save_path):
    """replots an H&E-stained image on a black background through simple thresholding rules

    args:
        image_path (str): file path that points to the underlying location of the TIF
        save_path (str): directory where output (replotted H&E image) will be written
    """

    image = cv2.imread(image_path)

    hsv_image = cv2.cvtColor(image, cv2.COLOR_BGR2HSV) # hue, saturation, value

    lower_purple = np.array([10, 13, 100])
    upper_purple = np.array([160, 255, 255])

    mask = cv2.inRange(hsv_image, lower_purple, upper_purple)

    purple_regions = cv2.bitwise_and(image, image, mask = mask)

    tiff.imwrite(os.path.join(save_path, os.path.basename(image_path)), cv2.cvtColor(purple_regions, cv2.COLOR_BGR2RGB))