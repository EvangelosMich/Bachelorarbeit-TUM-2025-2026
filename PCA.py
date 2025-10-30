import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from fancyimpute import IterativeSVD
import numpy as np

def nomralPCA():
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


def PCAISVD():
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

def PCAOALS():
    df = pd.read_csv("Dataset/S3(C)Expression.csv")
    df = df.set_index("ID").T
    D = df.to_numpy()
    R,C = D.shape
    N = 2
    PMatrix = np.random.random((N,C))
    TMatrix = np.zeros((R,N))
    tol = 1e-6
    max_iter = 200
    prev_error = np.inf

    for iteration in range(max_iter):
        for i in range(R):
            row = D[i,:]
            mask = ~np.isnan(row)
            if np.any(mask):
                P_sub = PMatrix[:,mask]
                d_sub = row[mask]
                
                TMatrix[i, :] = np.linalg.lstsq(P_sub.T,d_sub,rcond=None)[0]

        Q, _ = np.linalg.qr(TMatrix)
        TMatrix = Q

        for j in range(C):
            col = D[:,j]
            mask = ~np.isnan(col)
            if np.any(mask):
                T_sub = TMatrix[mask,:]
                d_sub = col[mask]

                PMatrix[:, j] = np.linalg.lstsq(T_sub, d_sub,rcond=None)[0]

        Qp,_ = np.linalg.qr(PMatrix.T)
        PMatrix = Qp.T        

        PMatrix = PMatrix/ np.linalg.norm(PMatrix,axis=1,keepdims=True)


        diff_sum = 0.0
        count = 0

        for i in range(R):
            for j in range(C):
                if not np.isnan(D[i,j]):
                    diff_sum += (D[i,j]-TMatrix[i,:] @ PMatrix[:,j]) ** 2
                    count+=1
        error = np.sqrt(diff_sum/count) if count > 0 else 0


        if abs(prev_error - error) < tol:
            print(f"Converged at iteration {iteration+1} with error {error:.6f}") 
            break

        prev_error = error


    else:
        print(f"Reached max iteration {max_iter}")           



    return TMatrix,PMatrix,df.index
    



def main():
    Tals,Pals,target = PCAOALS()

    reduced_df = pd.DataFrame(Tals,columns = ["PC1","PC2"])
    reduced_df["target"] = target.values

    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"],reduced_df["PC2"],s=50)

    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt,(reduced_df["PC1"][i],reduced_df["PC2"][i]),fontsize = 8)


    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("OALS PCA")
    plt.show()

    nomralPCA()

if __name__ == "__main__":
    main()    