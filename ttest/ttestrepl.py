import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv("Features_Matrix_Python.csv")



plt.figure(figsize=(8,6))
plt.scatter(-df["ttest.rapa_vs_VC.meandiff"], df["ttest.VC_vs_Rapa.logfdr"], s=50)

plt.xlabel("-MeanDiff")
plt.ylabel("Log.fdr")
plt.title("Comparison Ttest meandiff vs -fdr")
plt.show()