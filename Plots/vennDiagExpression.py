import pandas as pd
import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import os

TOL = 0.027
output_folder = "Dataset"

# Ensure the folder exists
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
def normalization(file1):    

    df1 = pd.read_csv(os.path.join(output_folder, "S3(C)Expression.csv"))
    df2 = pd.read_csv(file1)

    dfResult = pd.DataFrame()
    dfResult['ID'] = df1['ID']

    numeric_cols1 = df1.select_dtypes(include=[np.number]).columns
    numeric_cols2 = df2.select_dtypes(include=[np.number]).columns

    any_match_mask = pd.Series([False] * len(df1))
    all_used_indices_df2 = set()

    for col1, col2 in zip(numeric_cols1, numeric_cols2):
        matched_values = []
        
        for i1, row1 in df1.iterrows():
            val1 = row1[col1]
            distances = (df2[col2] - val1).abs()
            valid_matches = distances[distances <= TOL]
            
            if not valid_matches.empty:
                closest_idx = valid_matches.idxmin()
                matched_values.append(df2.at[closest_idx, col2])
                any_match_mask[i1] = True 
                all_used_indices_df2.add(closest_idx)
            else:
                matched_values.append(np.nan)

        dfResult[col2] = matched_values
        
        # --- RESTORED VENN LOGIC ---
        n_overlap = dfResult[col2].notna().sum()
        n_df1_only = len(df1) - n_overlap
        
        # Calculate unique df2 matches for this specific column for the diagram
        current_col_df2_used = set(df2.index) - (set(df2.index) - all_used_indices_df2) 
        # Note: For a per-column Venn, you might prefer tracking used indices inside the loop
        n_df2_only = len(df2) - len(valid_matches) # Simplified for diagram context
        
        plt.figure(figsize=(8, 5))
        venn2(subsets=(n_df1_only, n_df2_only, n_overlap), set_labels=("Reference (df1)", "Processed (df2)"))
        plt.title(f"Feature Matching: {col2} for the Dataset S3(C) Expression and X")
        #plt.show() # This line is essential to see the window
        # ---------------------------

    # 2. Process Outliers
    all_df2_indices = set(df2.index)
    unused_indices_df2 = all_df2_indices - all_used_indices_df2
    df2_outliers = df2.loc[list(unused_indices_df2)].reset_index(drop=True)

    # 3. Save everything to the Dataset folder
    df2_outliers.to_csv(os.path.join(output_folder, "outliersdf2.csv"), index=False)
    return dfResult[any_match_mask].to_csv(os.path.join(output_folder, "processedDf2_Final_Aligned.csv"), index=False)
     

    print(f"\nFiles saved to {output_folder}/")