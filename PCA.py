import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from fancyimpute import IterativeSVD
import numpy as np

def normalPCA():
    df = pd.read_csv("Dataset/S3(C)Expression.csv")
    df = df.set_index("ID").T
    # total_nan = df.isnull().sum().sum()
    # total_values = df.size
    # overall_percentage = (total_nan / total_values) * 100
    # print(f"Overall percentage of null values: {overall_percentage:.2f}%")
    target = df.index
    df = df.loc[:, df.isnull().mean() < 0.5]        #filter out rows and columns that have a NaN 
    df = df.loc[ df.isnull().mean(axis=1) < 0.5,:]  #ratio of bigger than .5 (can be changed)

    df = df.fillna(df.median(axis=0))               #Fill rest with median values of columns
    X_std = StandardScaler().fit_transform(df)      #Standarize is not used

    pca = PCA(n_components=2)                   # PCA with 2 components
    vecs = pca.fit_transform(X_std)                 


    reduced_df = pd.DataFrame(vecs,columns=["PC1","PC2"])
    reduced_df["PC1"] = -reduced_df["PC1"] #negate PC1 to get a more accurate model

    reduced_df["target"] = target.values




    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


    for i, txt in enumerate(reduced_df["target"]): #Name the values in order to have the clusters
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.xlabel(f'PC1({pca.explained_variance_ratio_[0]})')
    plt.ylabel(f'PC2({pca.explained_variance_ratio_[1]})')
    plt.show()


def PCAISVD():
    df = pd.read_csv("Dataset/S3(C)Expression.csv")
    df = df.set_index("ID").T
    #df = df.fillna(df.mean(),axis=0) 
    target = df.index

    imputer = IterativeSVD(rank=4) #Library by fancyimpute
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
    # Load the dataset and prepare it
    df = pd.read_csv("Dataset/S3(C)Expression.csv")
    df = df.set_index("ID").T           # transpose so rows = samples, columns = variables
    D = df.to_numpy()
    R, C = D.shape                      # R = number of samples, C = number of features
    N = 2                               # number of principal components to compute

    # Initialize P (loadings) randomly and T (scores) as zeros
    PMatrix = np.random.random((N, C))
    TMatrix = np.zeros((R, N))

    # Convergence settings
    tol = 1e-6
    max_iter = 200
    prev_error = np.inf

    # Main O-ALS loop
    for iteration in range(max_iter):

        # Step 1: Update T (scores) row by row using least squares
        for i in range(R):
            row = D[i, :]
            mask = ~np.isnan(row)  # only use available values (ignore NaNs)
            if np.any(mask):
                P_sub = PMatrix[:, mask]
                d_sub = row[mask]
                # solve for t(i,:) using least squares
                TMatrix[i, :] = np.linalg.lstsq(P_sub.T, d_sub, rcond=None)[0]

        # Orthogonalize T using QR decomposition (Gram-Schmidt)
        Q, _ = np.linalg.qr(TMatrix)
        TMatrix = Q

        # Step 2: Update P (loadings) column by column using least squares
        for j in range(C):
            col = D[:, j]
            mask = ~np.isnan(col)
            if np.any(mask):
                T_sub = TMatrix[mask, :]
                d_sub = col[mask]
                # solve for p(:,j)
                PMatrix[:, j] = np.linalg.lstsq(T_sub, d_sub, rcond=None)[0]

        # Orthogonalize P as well
        Qp, _ = np.linalg.qr(PMatrix.T)
        PMatrix = Qp.T

        # Normalize P rows so each loading vector has unit length
        PMatrix = PMatrix / np.linalg.norm(PMatrix, axis=1, keepdims=True)

        # Compute reconstruction error (only over observed entries)
        diff_sum = 0.0
        count = 0
        for i in range(R):
            for j in range(C):
                if not np.isnan(D[i, j]):
                    diff_sum += (D[i, j] - TMatrix[i, :] @ PMatrix[:, j]) ** 2
                    count += 1
        error = np.sqrt(diff_sum / count) if count > 0 else 0

        # Check for convergence (small change in error)
        if abs(prev_error - error) < tol:
            print(f"Converged at iteration {iteration+1} with error {error:.6f}")
            break

        prev_error = error

    else:
        print(f"Reached max iteration {max_iter}")

    # Return the scores, loadings, and sample labels
    return TMatrix, PMatrix, df.index
    



def main():
    T, P, target = PCAOALS()
    
    # Create DataFrame similar to your normalPCA output
    reduced_df = pd.DataFrame(P, columns=["PC1", "PC2"])
    reduced_df["PC1"] = -reduced_df["PC1"]  # same sign adjustment
    reduced_df["target"] = target.values

    # Plot
    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)

    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.xlabel("O-ALS PC1")
    plt.ylabel("O-ALS PC2")
    plt.title("PCA via Orthogonalized ALS (missing-value aware)")
    plt.show()

if __name__ == "__main__":
    main()    