import numpy as np
import pandas as pd
from fancyimpute import IterativeSVD
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


df = pd.read_csv("Dataset/S3(C)Expression.csv")
df = df.set_index("ID").T
df = df.fillna(df.mean(),axis=0)
target = df.index

imputer = IterativeSVD(rank = 5)
X_completed = imputer.fit_transform(df.values)

pca = PCA(n_components=2)
vecs = pca.fit_transform(X_completed)


reduced_df = pd.DataFrame(vecs,columns=["PC1","PC2"])
reduced_df["PC1"] = -reduced_df["PC1"]

reduced_df["target"] = target.values




plt.figure(figsize=(8,6))
plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


for i, txt in enumerate(reduced_df["target"]):
    plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

plt.xlabel(f'PC1({pca.explained_variance_ratio_[0]})')
plt.ylabel(f'PC2({pca.explained_variance_ratio_[1]})')
plt.show()






