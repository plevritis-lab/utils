import argparse
import logging
import numpy as np
import os
from skimage.segmentation import find_boundaries
from tifffile import imread, imwrite
from utils import extract_proteomic_panel
import warnings
import xml.etree.ElementTree as ET

warnings.filterwarnings("ignore", category = FutureWarning, message = ".*torch.load.*weights_only=False.*") # disable cellpose warning

# disable mesmer warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.getLogger().handlers = []
logging.basicConfig(level = logging.ERROR)

def construct_pseudochannel(image, segment_channel):
    """constructs a pseudochannel that merges intensities from a user-provided list of markers

    args:
        image (array): loaded matrix representation of the image of shape (c, y, x); \
                       introduces the pseudochannel at the end of the channels dimension, resulting in a new image of shape (c + 1, y, x)
        segment_channel (list): list of ordinal representations of the segment channels (0-indexed)
    """

    pseudochannel = np.sum(image[segment_channel, :, :], axis = 0)
    pseudochannel = np.clip(pseudochannel, 0, np.iinfo(image.dtype).max).astype(image.dtype)
    pseudochannel = pseudochannel[np.newaxis, :, :]

    return np.concatenate((image, pseudochannel), axis = 0)

def randomly_sample_patches(image, nuclear_channel, segment_channel, patch_size = (1000, 1000), 
                                num_patches = 5, intensity_threshold = 60, patch_threshold = 0.4):
    """selects and segments certain regions of the image with enough marker coverage for performance testing

    args:
        image (array): loaded matrix representation of the image of shape (c, y, x)
        nuclear_channel (int): ordinal representation of the nuclear channel (0-indexed) 
        segment_channel (int): ordinal representation of the segment channel (0-indexed)
        patch_size (tuple, optional): region size of the image to process of shape (y, x); \
                                      defaults to (1000, 1000)
        num_patches (int, optional): number of image regions to process; \
                                     defaults to 5
        intensity_threshold (int, optional): percentage that a particular patch pixel must exceed (in relation to the entire image intensity distribution) to contribute to criterion ranking; \
                                             defaults to 60
        patch_threshold (float, optional): percentage that a particular image patch must exceed to be considered for segmentation; \
                                           in particular, the sum of pixels with fluorescent intensity exceeding the intensity threshold in both the nuclear and segment channels in a patch \
                                               must exceed patch_size[0] * patch_size[1] * patch_threshold; \
                                           defaults to 0.4
    """

    y, x = image.shape[1], image.shape[2]

    upper_y = y - patch_size[0]
    upper_x = x - patch_size[1]

    patches = []

    nuclear_intensity_criterion = np.percentile(image[nuclear_channel, :, :], intensity_threshold)
    segment_intensity_criterion = np.percentile(image[segment_channel, :, :], intensity_threshold)

    while len(patches) < num_patches:
        sampled_y = np.random.choice(upper_y)
        sampled_x = np.random.choice(upper_x)

        patch = image[:, sampled_y : sampled_y + patch_size[0],
                         sampled_x : sampled_x + patch_size[1]]
        
        nuclear_channel_intensities = patch[nuclear_channel, :, :]
        segment_channel_intensities = patch[segment_channel, :, :]

        nuclear_criterion_count = np.count_nonzero(nuclear_channel_intensities >= nuclear_intensity_criterion)
        segment_criterion_count = np.count_nonzero(segment_channel_intensities >= segment_intensity_criterion)

        nuclear_criteria = nuclear_criterion_count / (patch_size[0] * patch_size[1])
        segment_criteria = segment_criterion_count / (patch_size[0] * patch_size[1])

        if nuclear_criteria > patch_threshold and segment_criteria > patch_threshold:
            patches.append(patch)

    return patches

def compress_channels(images, nuclear_channel, segment_channel, save_path):
    """compresses the multichannel QPTIFF into an RGB image, with G encapsulating the segment channel and B encapsulating the nuclear channel

    args:
        images (list): list of loaded images of shape (c, y, x)
        nuclear_channel (int): ordinal representation of the nuclear channel (0-indexed)
        segment_channel (tuple): tuple of form (channel name, ordinal representation of the segment channel (0-indexed))
        save_path (str): parent destination where output (RGB images) will be written under their own subdirectory (original/{segment_channel[0]})
    """

    compressed_images = [None] * len(images)

    segment_channel_name, segment_channel = segment_channel
    
    save_path = os.path.join(save_path, "original")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok = True)

    for i, image in enumerate(images):
        compressed_image = np.zeros((3, image.shape[1], image.shape[2]), dtype = image.dtype)
        compressed_image[1, :, :] = image[segment_channel, :, :] # color the segment channel G
        compressed_image[2, :, :] = image[nuclear_channel, :, :] # color the nuclear channel B

        imwrite(os.path.join(save_path, f"image_{i + 1}.tif"), np.transpose(compressed_image, (1, 2, 0)))
        compressed_images[i] = compressed_image
    
    return compressed_images

def apply_cellpose(compressed_images, segment_channel_name, save_path):
    """applies a pretrained cellpose model (cyto3) to segment images of shape (y, x, c) or (c, y, x)

    args:
        compressed_images (list): list of RGB images of shape (3, y, x), where G is the segment channel and B is the nuclear channel
        segment_channel_name (str): name of the protein used as the segment channel
        save_path (str): parent destination where output (RGB images, segmentation masks) will be written under their own subdirectory (cellpose/segment_channel_name)
    """

    from cellpose import io, models, utils

    model = models.Cellpose(model_type = "cyto3")

    save_path = os.path.join(save_path, "cellpose")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok = True)

    for i, compressed_image in enumerate(compressed_images):
        masks, flows, _, diams = model.eval(compressed_image, [2, 3], diameter = None)
        compressed_image[0, :, :] = utils.masks_to_outlines(masks) * np.iinfo(compressed_image.dtype).max # color the segmentation mask R

        imwrite(os.path.join(save_path, f"image_{i + 1}.tif"), np.transpose(compressed_image, (1, 2, 0)))
        io.masks_flows_to_seg(compressed_image, masks, flows, os.path.join(save_path, f"image_{i + 1}"), diams, [2, 3])

def apply_mesmer(compressed_images, segment_channel_name, save_path):
    """applies a pretrained deepcell model (mesmer) to segment images of shape (batch, x, y, c)

    args:
        compressed_images (list): list of RGB images of shape (3, y, x), where G is the segment channel and B is the nuclear channel
        segment_channel_name (str): name of the protein used as the segment channel
        save_path (str): parent destination where output (RGB images, segmentation masks) will be written under their own subdirectory (mesmer/segment_channel_name)
    """

    from deepcell.applications import Mesmer
    from deepcell.utils.plot_utils import make_outline_overlay, create_rgb_image

    os.environ.update({"DEEPCELL_ACCESS_TOKEN" : "HcI12JDz.jBnsD080yqdai9s5f0LHZodVdEUpPimh"})

    model = Mesmer()

    save_path = os.path.join(save_path, "mesmer")
    save_path = os.path.join(save_path, segment_channel_name)
    os.makedirs(save_path, exist_ok = True)

    for i, compressed_image in enumerate(compressed_images):
        reshaped_image = np.transpose(compressed_image, (2, 1, 0)).copy()
        
        reshaped_image = reshaped_image[:, :, 1:] # remove empty R channel
        reshaped_image = reshaped_image[..., ::-1] # swap G (segment channel) and B (nuclear channel)

        reshaped_image = np.expand_dims(reshaped_image, axis = 0)

        segmentation_predictions = model.predict(reshaped_image, batch_size = 1, image_mpp = 0.377, compartment = "whole-cell")[0, :, :, 0]
        segmentation_predictions = np.transpose(segmentation_predictions, (1, 0))

        reshaped_image = reshaped_image[0, :, :, :]
        reshaped_image = np.stack((np.zeros((reshaped_image.shape[0], reshaped_image.shape[1]), dtype = reshaped_image.dtype), 
                                       reshaped_image[:, :, 1], reshaped_image[:, :, 0]), axis = -1)
        
        reshaped_image = np.transpose(reshaped_image, (1, 0, 2))

        located_boundaries = find_boundaries(segmentation_predictions, connectivity = 1, mode = "inner")
        reshaped_image[:, :, 0][located_boundaries > 0] = np.iinfo(reshaped_image.dtype).max

        imwrite(os.path.join(save_path, f"image_{i + 1}.tif"), reshaped_image)
        np.save(os.path.join(save_path, f"image_{i + 1}_seg.npy"), segmentation_predictions)

def overlay_masks(compressed_images, cellpose_segment_channel_name, mesmer_segment_channel_name, display_segment_channel_name, load_path, save_path):
    """overlays segmentation masks generated by both cellpose and mesmer onto the same image

    args:
        compressed_images (list): list of RGB images of shape (3, y, x), where G is the display segment channel and B is the nuclear channel
        cellpose_segment_channel_name (str): segmentation channel name used by cellpose
        mesmer_segment_channel_name (str): segmentation channel name used by mesmer
        display_segment_channel_name (str): segmentation channel name used for visual display
        load_path (str): base folder where segmentation masks (.npy files) are stored; \
                         given the symmetry in how mesmer and cellpose outputs are saved, only specify the topmost output directory (e.g., patches or full)
        save_path (str): parent destination where output (overlaid image masks) will be written under their own subdirectory (combined_algorithm_overlays)
    """

    save_path = os.path.join(save_path, "combined_algorithm_overlays")
    os.makedirs(save_path, exist_ok = True)

    for i, compressed_image in enumerate(compressed_images):
        overlaid_image = compressed_image.copy()

        try:
            cellpose_segmentation = np.load(os.path.join(load_path, f"cellpose/{cellpose_segment_channel_name}/image_{i + 1}_seg.npy"), allow_pickle = True).item()["outlines"]
            mesmer_segmentation = np.load(os.path.join(load_path, f"mesmer/{mesmer_segment_channel_name}/image_{i + 1}_seg.npy"), allow_pickle = True)

        except FileNotFoundError:
            print("\nERROR: please double-check that your provided segmentation masks are in their specified locations\n"
                  f"currently using cellpose segmentation mask path: {os.path.join(load_path, f'cellpose/{cellpose_segment_channel_name}/image_{i + 1}_seg.npy')}\n"
                  f"currently using mesmer segmentation mask path: {os.path.join(load_path, f'mesmer/{mesmer_segment_channel_name}/image_{i + 1}_seg.npy')}")
            
            return

        mesmer_boundaries = find_boundaries(mesmer_segmentation, connectivity = 1, mode = "inner")

        yellow_mask = cellpose_segmentation > 0 
        red_mask = mesmer_boundaries > 0 

        # TODO - consider handcrafting a legend for the overlaid masks, with labels for the segment channel as well as red and yellow outlines (this requires a nontrivial fix without using matplotlib)
        
        overlaid_image[0, :, :][yellow_mask] = np.iinfo(overlaid_image.dtype).max
        overlaid_image[1, :, :][yellow_mask] = np.iinfo(overlaid_image.dtype).max # color cellpose mask Y 
        overlaid_image[0, :, :][red_mask] = np.iinfo(overlaid_image.dtype).max # color mesmer mask R
        overlaid_image[1, :, :][red_mask & yellow_mask] = int(np.iinfo(overlaid_image.dtype).max * 2 / 3) # color intersection O

        imwrite(os.path.join(save_path, f"image_{i + 1}_cellpose_{cellpose_segment_channel_name}_mesmer_{mesmer_segment_channel_name}_display_{display_segment_channel_name}.tif"),
                   np.transpose(overlaid_image, (1, 2, 0)))

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    # TODO - make this argument list a little nicer (have something like 'git commit <args>' and 'git merge <args>') but for this, i.e., 'apply_segmentation segment <args>' and 'apply_segmentation overlay <args>'

    parser = argparse.ArgumentParser(description = "interface for applying cellular segmentation algorithms, mesmer and cellpose; \
                                                    it is recommended to choose a cytoplasmic + nuclear combination for cellpose; \
                                                    it is recommended to choose a membrane + nuclear combination for mesmer")
        
    parser.add_argument("-i", "--image_path", help = "file path that points to the underlying location of the QPTIFF or TIF")
    parser.add_argument("-p", "--panel_path", help = "file path that points to the underlying location of the channel_names.txt file; \
                                                      this argument is optional and should only be used when metadata parsing fails")
    parser.add_argument("-n", "--nuclear_channel", help = "name of the nuclear channel to be used in segmentation OR a display channel if overlay_masks is toggled")
    parser.add_argument("-s", "--segment_channel", help = "name of the cytoplasm or membrane channel(s) to be used in segmentation OR a display channel if overlay_masks is toggled; \
                                                           if multiple channels are desired, please supply a comma separated list with no spaces, which will construct a pseudochannel that merges their intensities")

    parser.add_argument("--apply_cellpose", action = "store_true", help = "segments the image using cellpose")
    parser.add_argument("--apply_mesmer", action = "store_true", help = "segments the image using memser")
    parser.add_argument("--overlay_masks", action = "store_true", help = "overlays the segmentation masks produced by both cellpose and mesmer onto the original image for visual inspection")
    
    parser.add_argument("-c", "--cellpose_channel", help = "name of the protein channel used by cellpose to generate the image segmentation mask; \
                                                            only relevant if overlay_masks is toggled")
    parser.add_argument("-m", "--mesmer_channel", help = "name of the protein channel used by mesmer to generate the image segmentation mask; \
                                                          only relevant if overlay_masks is toggled")
    
    parser.add_argument("--debug", action = "store_true", help = "selects and segments certain regions of the image with enough marker coverage for performance testing")

    return parser.parse_args()
    
def main():
    arguments = parse_arguments()

    image_path = arguments.image_path
    panel_path = arguments.panel_path
    nuclear_channel = arguments.nuclear_channel
    segment_channel = arguments.segment_channel

    segment_channel = segment_channel.split(",")

    np.random.seed(42)

    print("\nprocessing sample", os.path.basename(os.path.splitext(image_path)[0]))

    try:
        panel = extract_proteomic_panel(image_path)
    
    except (AttributeError, KeyError, ET.ParseError):
        print("\nERROR: unable to extract marker metadata from the provided file; "
                  "please provide a channel_names.txt file instead")
        
        if panel_path:
            print(f"INFO: extracting protein panel from {panel_path}")
            panel = extract_proteomic_panel(image_path, panel_path)
        
        else:
            return

    try:
        nuclear_channel_index = panel.index(nuclear_channel)
        segment_channel_index = [panel.index(s) for s in segment_channel]
    
    except ValueError:
        print(f"\nERROR: please ensure that the provided nuclear ({nuclear_channel}) and segment " 
                  f"({segment_channel}) channel names are in your protein panel")
        return
    
    image = imread(image_path) # (c, y, x)

    if len(segment_channel) > 1:
        image = construct_pseudochannel(image, segment_channel_index)
        
        segment_channel = " + ".join(segment_channel)
        segment_channel_index = image.shape[0] - 1

    else:
        segment_channel = segment_channel[0]
        segment_channel_index = segment_channel_index[0]

    base_directory = os.path.dirname(os.path.dirname(image_path))

    if arguments.debug:
        save_path = os.path.join(base_directory, "patches")

        patches = randomly_sample_patches(image, nuclear_channel_index, segment_channel_index, intensity_threshold = 90, patch_threshold = 0.2)
        compressed_images = compress_channels(patches, nuclear_channel_index, (segment_channel, 
                                                                               segment_channel_index), save_path)

    else:
        save_path = os.path.join(base_directory, "full")

        compressed_images = compress_channels([image], nuclear_channel_index, (segment_channel, 
                                                                               segment_channel_index), save_path)
        
    if arguments.apply_cellpose:
        apply_cellpose(compressed_images, segment_channel, save_path)
    
    if arguments.apply_mesmer:
        apply_mesmer(compressed_images, segment_channel, save_path)
    
    if arguments.overlay_masks:
        cellpose_channel = arguments.cellpose_channel
        mesmer_channel = arguments.mesmer_channel

        overlay_masks(compressed_images, cellpose_channel, mesmer_channel, segment_channel, save_path, save_path)

main()