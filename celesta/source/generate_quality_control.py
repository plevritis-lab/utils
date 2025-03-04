import argparse
import os
import pandas as pd

def generate_quality_control_spreadsheet(clinical_data_path, panel_path, save_path):
    """generates a quality control spreadsheet for spatial proteomics data

    args:
        clinical_data_path (str): path to a CSV file containing clinical data
        panel_path (str): path to a TXT file containing the list of proteins in the panel
        save_path (str): path where the generated quality control spreadsheet will be saved
    """

    clinical_data = pd.read_csv(clinical_data_path)
    with open(panel_path, "r") as panel:
        protein_panel = [m.rstrip() for m in panel]

    grouped = clinical_data.groupby(["TMA", "TMA_PART"])

    output_data = []
    previous_microarray_id = None

    for (tma, part), group in grouped:
        microarray_id = f"{tma}_{part}"

        if previous_microarray_id and previous_microarray_id != microarray_id:
            output_data.append([""] * (2 + len(protein_panel)))

        for _, row in group.iterrows():
            sample_id = f"reg{int(row['CORE_IMAGE_ID']):03d}"
            output_data.append([microarray_id, sample_id] + [""] * len(protein_panel))

        previous_microarray_id = microarray_id

    columns = ["MICROARRAY", "SAMPLE"] + protein_panel

    quality_control_spreadsheet = pd.DataFrame(output_data, columns = columns)
    quality_control_spreadsheet.to_csv(save_path + "_quality_control.csv", index = False)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for generating quality control spreadsheets")

    parser.add_argument("-d", "--clinical_data", help = "path to a clinical data spreadsheet that contains columns 'TMA', 'TMA_PART', and 'CORE_IMAGE_ID'")
    parser.add_argument("-p", "--panel_path", help = "path that points to the underlying location of the channel_names.txt file")
    parser.add_argument("-s", "--save_path", help = "path where the generated quality control spreadsheet will be saved; \
                                                     this argument will be suffixed by _quality_control.csv")
    
    return parser.parse_args() 
    
def main():
    arguments = parse_arguments()

    clinical_data = arguments.clinical_data
    panel_path = arguments.panel_path
    save_path = arguments.save_path

    os.makedirs(os.path.dirname(save_path), exist_ok = True)

    generate_quality_control_spreadsheet(clinical_data, panel_path, save_path)

main()