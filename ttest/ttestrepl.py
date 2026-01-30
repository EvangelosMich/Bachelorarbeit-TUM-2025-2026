import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv("Dataset/Features_Aligned_Final.csv")



plt.figure(figsize=(8,6))
plt.scatter(-np.log10(df["ttest.Rapamycin_vs_Control.mean.diff"].values), df["ttest.Rapamycin_vs_Control.log.fdr"].values, s=50)



plt.xlabel("-MeanDiff")
plt.ylabel("Log.fdr")
plt.title("Comparison Ttest meandiff vs -fdr for the normalized features sheet")
plt.show()