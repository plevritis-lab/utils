import argparse
import os
import pandas as pd

def generate_initial_thresholds(clinical_data_path, signature_matrix_path, save_path):
    """generates thresholds spreadsheets for spatial proteomics data

    args:
        clinical_data_path (str): path to a CSV file containing clinical data
        signature_matrix_path (str): path to a CSV file containing a signature matrix
        save_path (str): path where the generated thresholds spreadsheets will be saved
    """

    clinical_data = pd.read_csv(clinical_data_path)
    signature_matrix = pd.read_csv(signature_matrix_path)

    for _, row in clinical_data.iterrows():
        core_image_id = f"reg{int(row['CORE_IMAGE_ID']):03d}"
        identifier = f"{row['TMA']}_{row['TMA_PART']}_{core_image_id}"
        
        thresholds = pd.DataFrame({
            "CELL_TYPE": signature_matrix["CELL_TYPE"],
            "ANCHOR": [0.7] * len(signature_matrix),
            "INDEX": [0.5] * len(signature_matrix)
        })
                
        thresholds.to_csv(os.path.join(save_path, f"{identifier}_thresholds.csv"), index = False)
    
def update_thresholds(thresholds_path, signature_matrix_path):
    """updates existing threshold files with missing cell types from the signature matrix
    
    agrs:
        thresholds_path (str): path to directory containing threshold CSV files
        signature_matrix_path (str): path to the signature matrix CSV file
    """

    signature_matrix = pd.read_csv(signature_matrix_path)

    for file in os.listdir(thresholds_path):
        if not file.endswith("_thresholds.csv"):
            continue
            
        thresholds = pd.read_csv(os.path.join(thresholds_path, file))
        
        rows = []
        for cell_type in signature_matrix["CELL_TYPE"]:
            if cell_type in thresholds["CELL_TYPE"].values:
                row = thresholds[thresholds["CELL_TYPE"] == cell_type].iloc[0]
                rows.append({
                    "CELL_TYPE": cell_type,
                    "ANCHOR": row["ANCHOR"],
                    "INDEX": row["INDEX"]
                })

            else:
                rows.append({
                    "CELL_TYPE": cell_type,
                    "ANCHOR": 0.7,
                    "INDEX": 0.5
                })
        
        updated_thresholds = pd.DataFrame(rows)
        updated_thresholds.to_csv(os.path.join(thresholds_path, file), index = False)

def parse_arguments():
    """parses several command line arguments provided by the user (use --help to see the full list)"""

    parser = argparse.ArgumentParser(description = "interface for generating thresholds spreadsheets")

    parser.add_argument("-d", "--clinical_data", help = "path to a clinical data spreadsheet that contains columns 'TMA', 'TMA_PART', and 'CORE_IMAGE_ID'")
    parser.add_argument("-m", "--signature_matrix", help = "path that points to the underlying location of the signature matrix")
    parser.add_argument("-s", "--save_path", help = "path where the generated thresholds spreadsheets will be saved")

    parser.add_argument("--update", action = "store_true", help = "updates existing threshold files with missing cell types from the signature matrix")
    
    return parser.parse_args() 
    
def main():
    arguments = parse_arguments()

    clinical_data = arguments.clinical_data
    signature_matrix = arguments.signature_matrix
    save_path = arguments.save_path

    os.makedirs(save_path, exist_ok = True)

    update_thresholds(save_path, signature_matrix) if arguments.update else \
        generate_initial_thresholds(clinical_data, signature_matrix, save_path)

main()