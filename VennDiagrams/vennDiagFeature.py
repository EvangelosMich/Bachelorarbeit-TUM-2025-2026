import pandas as pd
import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import os

MZTOL = 1.5
RTTOL = 5.0
output_folder = "Dataset"

# 1. Load Data
df1 = pd.read_csv("Dataset/processedFeatureInfo.csv")
df2 = pd.read_csv("Dataset/S3(B)Feature.csv")

# Initialize the result dataframe with metadata from df1
dfResult = df1.copy()

# Tracking sets
matched_indices_df1 = []
used_indices_df2 = set()
df2_aligned_rows = []

# 2. Matching Logic
for i1, row1 in df1.iterrows():
    mz1 = row1["General.All.mzmed"]
    rt1 = row1["General.All.rtmed"]

    # Filter df2 for candidates within both tolerances
    mask = (
        (abs(df2["General|All|mzmed"] - mz1) <= MZTOL) &
        (abs(df2["General|All|rtmed"] - rt1) <= RTTOL)
    )
    matches = df2[mask]

    if not matches.empty:
        # If multiple matches, find the one with the smallest combined Euclidean distance
        # We normalize slightly or just use raw MZ difference if RT is stable
        distances = np.sqrt((matches["General|All|mzmed"] - mz1)**2 + (matches["General|All|rtmed"] - rt1)**2)
        closest_idx = distances.idxmin()
        
        # Store index for Venn and data for merging
        matched_indices_df1.append(i1)
        used_indices_df2.add(closest_idx)
        
        # Grab the full row from df2 to append to df1's metadata later
        df2_aligned_rows.append(df2.loc[closest_idx])
    else:
        # No match found for this df1 feature
        df2_aligned_rows.append(pd.Series([np.nan] * len(df2.columns), index=df2.columns))

# 3. Construct Final Aligned DataFrame
# We join the columns of df1 with the matched columns of df2
df2_matched_df = pd.DataFrame(df2_aligned_rows).reset_index(drop=True)
# To avoid column name collisions, you might want to suffix df2 columns
dfFinal = pd.concat([df1, df2_matched_df.add_suffix('_df2')], axis=1)

# 4. Filter for only the successful matches (like any_match_mask)
dfFinal_Filtered = dfFinal.iloc[matched_indices_df1]

# 5. Save Files
dfFinal_Filtered.to_csv(os.path.join(output_folder, "Features_Aligned_Final.csv"), index=False)

# Outliers (Features in df2 that never matched a reference in df1)
unused_indices_df2 = set(df2.index) - used_indices_df2
df2_outliers = df2.loc[list(unused_indices_df2)]
df2_outliers.to_csv(os.path.join(output_folder, "outliers_features_df2.csv"), index=False)

# 6. Venn Diagram
n_overlap = len(matched_indices_df1)
n_df1_only = len(df1) - n_overlap
n_df2_only = len(df2) - len(used_indices_df2)

plt.figure(figsize=(8, 5))
venn2(subsets=(n_df1_only, n_df2_only, n_overlap), set_labels=("Ref Features (df1)", "Target Features (df2)"))
plt.title("MZ/RT Feature Alignment Results")
plt.show()

print(f"Alignment complete. Saved matched features and {len(df2_outliers)} outliers to {output_folder}/")


