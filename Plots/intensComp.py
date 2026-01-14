import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




df1 = pd.read_csv("Dataset/processedDf2_Final_Aligned.csv")
df2 = pd.read_csv("Dataset/outliersdf2.csv")
df3 = pd.read_csv("Dataset/ExpressionProcessed.csv")
df4 = pd.read_csv("Dataset/S3(C)Expression.csv")

numCol1 = df1.select_dtypes(include=[np.number]).columns
numCol2 = df2.select_dtypes(include=[np.number]).columns
numCol3 = df3.select_dtypes(include=[np.number]).columns
numCol4 = df4.select_dtypes(include=[np.number]).columns


for col1,col2,col3,col4 in zip(numCol1,numCol2,numCol3,numCol4):
    meanVal1 = df1[col1].mean()
    meanVal2 = df2[col2].mean()
    meanVal3 = df3[col3].mean()
    meanVal4 = df4[col4].mean()

    name = ['Processed-noOutlier','Outlier','newData','origiData']
    values = [meanVal1,meanVal2,meanVal3,meanVal4]

    bars = plt.bar(name,values)
    plt.bar_label(bars,padding=3)
    plt.title(f"Comparison of average values of each column for{col1}")
    plt.show()