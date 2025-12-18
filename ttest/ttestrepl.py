import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv("Dataset/S3(B)Feature.csv")



plt.figure(figsize=(8,6))
plt.scatter(-df["ttest|Control_vs_Rapamycin|mean.diff"], df["ttest|Control_vs_Rapamycin|log.fdr"], s=50)

plt.xlabel("-MeanDiff")
plt.ylabel("Log.fdr")
plt.title("Comparison Ttest meandiff vs -fdr")
plt.show()