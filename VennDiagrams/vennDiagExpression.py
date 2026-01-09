import pandas as pd
import numpy as np
from matplotlib_venn import venn2,venn3
from matplotlib import pyplot as plt
from thefuzz import fuzz,process

TOL = 0.027



df1 = pd.read_csv("Dataset/S3(C)Expression.csv")
df2 = pd.read_csv("Dataset/ExpressionProcessed.csv")

df3 = []

numeric_cols1 = df1.select_dtypes(include=[np.number]).columns
numeric_cols2 = df2.select_dtypes(include=[np.number]).columns

for col1,col2 in zip(numeric_cols1,numeric_cols2):
    df1_matched_features = set()
    df2_matched_features = set()
    df2_indices_used = set()

    for i1,row1 in df1.iterrows():
        mz1 = row1[col1]
        

        matches = df2[
            (abs(df2[col2]-mz1)<=TOL) 
        ]

        
        if len(matches) > 0:
            df1_matched_features.add(i1)
            df2_indices_used.update(matches.index)
        


    df1_indices_used = set()    
    for i2,row2 in df2.iterrows():
        mz2 = row2[col2]

        matches = df1[(abs(df1[col1]-mz2)<=TOL)]

        if len(matches) > 0:
            df2_matched_features.add(i2)
            df1_indices_used.update(matches.index)        



    alldf2Indices = set(df2.index)
    alldf1Indices = set(df1.index)

    matching = df1.loc[list(df1_matched_features)]
    

    


    outlier_indicesdf2 = alldf2Indices - df2_indices_used
    outlier_indiceddf1 = alldf1Indices - df1_indices_used

    outlierdf1 = df1.loc[list(outlier_indiceddf1)]
    outliersdf2 = df2.loc[list(outlier_indicesdf2)]

    

    outlierdf1[col1] = outlierdf1[col1].fillna(0)

    mean_value_outlier_df1 = outlierdf1[col1].mean()
    mean_value_outlier_df2 = outliersdf2[col2].mean()
    mean_value_matches = matching[col1].mean()

    print(col1)
    print(mean_value_outlier_df1)
    print(mean_value_outlier_df2)
    print(mean_value_matches)


    n_overlap = len(df1_matched_features)
    n_df1_only = len(df1) - n_overlap
    n_df2_only = len(df2)-len(df2_indices_used)

    outlierdf1.to_csv("outliersDf1.csv",index=False)
    outliersdf2.to_csv("outliersDf2.csv",index=False)


    venn2(subsets=(n_df1_only,n_df2_only,n_overlap),set_labels=("XY"))
    plt.show()



