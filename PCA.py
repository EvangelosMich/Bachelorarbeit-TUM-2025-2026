import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler



df = pd.read_csv("Dataset/S3(C)Expression.csv")
df = df.set_index("ID").T
total_nan = df.isnull().sum().sum()
total_values = df.size
overall_percentage = (total_nan / total_values) * 100
print(f"Overall percentage of null values: {overall_percentage:.2f}%")
target = df.index
df = df.fillna(df.median(axis=0))

X_std = StandardScaler().fit_transform(df)

pca = PCA(n_components=2)
vecs = pca.fit_transform(df)


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




