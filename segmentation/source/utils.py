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
        with TiffFile(image_path) as tif:
            for page in tif.series[0].pages:
                image_metadata = page.tags["ImageDescription"].value
                image_metadata = ElementTree.fromstring(image_metadata)

                protein = image_metadata.find("Biomarker").text

                protein_panel.append(protein)

    return protein_panel