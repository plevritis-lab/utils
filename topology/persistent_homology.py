import pandas as pd
import os
import glob
import gudhi as gd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import re

assignments_directory = "$HOME/Downloads/example"
assignments_directory = os.path.expandvars(assignments_directory)

csv_files = glob.glob(os.path.join(assignments_directory, "*.csv"))

for csv_file in csv_files:
    filename = os.path.basename(csv_file)
    df = pd.read_csv(csv_file)

    # Extract the pattern from the filename
    pattern = re.search(r'(\d+_S\d+_reg\d+)', filename)
    if pattern:
        sample_id = pattern.group(1)
        print(f"Extracted sample ID: {sample_id}")
    else:
        sample_id = os.path.splitext(filename)[0]  # Fallback to filename without extension
        print(f"Could not extract sample ID pattern, using: {sample_id}")

    if sample_id == "521_S1_reg024":
        # Group by cell type
        cell_types = df['FINAL_CELL_TYPE'].unique()

        # Create a directory for saving the diagrams if it doesn't exist
        output_dir = os.path.expandvars("$HOME/Downloads/topology")
        os.makedirs(output_dir, exist_ok=True)

        # Process each cell type
        for cell_type in cell_types:
            print(f"Processing {cell_type} from {filename}")
            
            # Filter data for this cell type
            cell_df = df[df['FINAL_CELL_TYPE'] == cell_type]
            
            # Skip if not enough points
            if len(cell_df) < 3:
                print(f"Skipping {cell_type}, not enough points")
                continue
            
            # Extract coordinates for persistent homology
            points = cell_df[['X', 'Y']].values
            
            # Create a Rips complex
            rips_complex = gd.RipsComplex(points=points)
            simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
            
            # Compute persistence diagrams
            persistence = simplex_tree.persistence()
            
            # Plot the persistence diagram
            plt.figure(figsize=(10, 8))
            gd.plot_persistence_diagram(persistence)
            plt.title(f"{cell_type}")
            
            # Save the figure
            output_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_{cell_type}_ph_diagram.pdf")
            plt.savefig(output_path, dpi=300)
            # plt.savefig(output_path)
            plt.close()

            # Calculate persistence entropy
            def persistent_entropy(persistence):
                # Filter out features with infinite death time
                finite_pairs = [(birth, death) for dim, (birth, death) in persistence if death != float('inf')]
                
                if not finite_pairs:
                    return 0.0
                
                # Calculate lifetimes (death - birth)
                lifetimes = [death - birth for birth, death in finite_pairs]
                total_lifetime = sum(lifetimes)
                
                # Calculate normalized lifetimes (probability distribution)
                if total_lifetime == 0:
                    return 0.0
                
                p = [lt / total_lifetime for lt in lifetimes]
                
                # Calculate entropy
                return -sum(pi * np.log(pi) for pi in p if pi > 0)

            # Extract features from persistence diagram
            def extract_features(persistence):
                # Get features by dimension
                features_by_dim = defaultdict(list)
                for dim, (birth, death) in persistence:
                    if death != float('inf'):  # Skip infinite features
                        features_by_dim[dim].append(death - birth)  # Store lifetime
                
                # Calculate statistics for each dimension
                feature_vector = []
                
                # For dimensions 0 and 1
                for dim in [0, 1]:
                    lifetimes = features_by_dim[dim]
                    if lifetimes:
                        feature_vector.extend([
                            len(lifetimes),            # Number of features
                            np.mean(lifetimes),        # Mean lifetime
                            np.std(lifetimes) if len(lifetimes) > 1 else 0,  # Standard deviation
                            np.max(lifetimes) if lifetimes else 0,  # Max lifetime
                            np.sum(lifetimes)          # Total persistence
                        ])
                    else:
                        feature_vector.extend([0, 0, 0, 0, 0])  # Default values if no features
                
                # Add persistence entropy
                feature_vector.append(persistent_entropy(persistence))
                
                return feature_vector

            # Calculate and store feature vectors
            feature_vector = extract_features(persistence)
            print(f"Persistence feature vector for {cell_type}:")
            print(f"[Number of H0, Mean H0 lifetime, Std H0 lifetime, Max H0 lifetime, Total H0 persistence, "
                  f"Number of H1, Mean H1 lifetime, Std H1 lifetime, Max H1 lifetime, Total H1 persistence, "
                  f"Persistence Entropy]")
            print(feature_vector)

            # Define the output file path
            features_output_file = os.path.join(output_dir, f"{sample_id}_persistence_features.txt")

            # Create a header if the file doesn't exist
            if not os.path.exists(features_output_file) or os.path.getsize(features_output_file) == 0:
                with open(features_output_file, 'w') as f:
                    f.write("Sample_ID,Cell_Type,Num_H0,Mean_H0_Lifetime,Std_H0_Lifetime,Max_H0_Lifetime,Total_H0_Persistence,"
                            "Num_H1,Mean_H1_Lifetime,Std_H1_Lifetime,Max_H1_Lifetime,Total_H1_Persistence,Persistence_Entropy\n")

            # Append the feature vector with the cell type
            with open(features_output_file, 'a') as f:
                feature_str = ','.join([str(val) for val in feature_vector])
                f.write(f"{sample_id},{cell_type},{feature_str}\n")

            print(f"Feature vector for {cell_type} saved to {features_output_file}")

        break