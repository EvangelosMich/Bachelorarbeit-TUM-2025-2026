import pandas as pd
import numpy as np
from matplotlib_venn import venn2,venn3
from matplotlib import pyplot as plt
from thefuzz import fuzz,process

MZTOL = 1.5
RTTOL = 5.0



df1 = pd.read_csv("Dataset/S3(B)Feature.csv")
df2 = pd.read_csv("Dataset/processedFeatureInfo.csv")



matched_features = set()
df2_indices_used = set()

for i1,row1 in df1.iterrows():
    mz1 = row1["General|All|mzmed"]
    rt1 = row1["General|All|rtmed"]


    matches = df2[
        (abs(df2["General.All.mzmed"]-mz1)<=MZTOL) &
        (abs(df2["General.All.rtmed"]-rt1)<=RTTOL)
    ]


    if len(matches) > 0:
        matched_features.add(i1)
        df2_indices_used.update(matches.index)

n_overlap = len(matched_features)
n_df1_only = len(df1) - n_overlap
n_df2_only = len(df2)-len(df2_indices_used)

venn2(subsets=(n_df1_only,n_df2_only,n_overlap),set_labels=("File 1 Features,File2Features"))
plt.show()


