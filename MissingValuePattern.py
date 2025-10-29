import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

df = pd.read_csv("Dataset/S3(C)Expression.csv")

mask = df.isna()

col_na = mask.mean()
row_na = mask.mean(axis=1)


print(col_na.std(),row_na.std())

plt.imshow(mask,aspect="auto",cmap="gray")
plt.xlabel("Variables")
plt.ylabel("Samples")
#lt.show()

mask = df.isna().astype(int)
clusters = KMeans(n_clusters=3, random_state=0).fit_predict(mask)
print(pd.Series(clusters).value_counts())