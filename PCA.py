import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler



df = pd.read_csv("Dataset/S3(C)Expression.csv")
df = df.set_index("ID").T
target = df.index
df = df.fillna(df.mean(axis=0))
X_std = StandardScaler().fit_transform(df)

pca = PCA(n_components=2)
vecs = pca.fit_transform(df)


reduced_df = pd.DataFrame(vecs,columns=["PC1","PC2"])
reduced_df["target"] = target.values
print(pca.explained_variance_ratio_)


plt.figure(figsize=(8,6))
plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


for i, txt in enumerate(reduced_df["target"]):
    plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

plt.xlabel('PC1(42.59%)')
plt.ylabel('PC2(16.38%)')
plt.show()




