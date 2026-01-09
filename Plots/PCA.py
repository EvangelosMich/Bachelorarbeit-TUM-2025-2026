import pandas as pd
from sklearn.decomposition import PCA   
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from fancyimpute import IterativeSVD
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def normalPCA():
    # Load data and transpose so rows=samples, cols=features
    df = pd.read_csv("Dataset/S3(C)Expression.csv").set_index("ID").T

    # Log10 transform (keep values from exploding, avoid log(0) with +1)
    df = df.apply(lambda x: np.log10(x + 1) if np.issubdtype(x.dtype, np.number) else x)

    # Sample labels for plotting
    target = df.index

    # Fill missing values + remove zero-variance features
    df = df.fillna(0)
    df = df.loc[:, df.var(ddof=1) > 0]

    # Run PCA
    X = df.to_numpy(dtype=float)
    pca = PCA(n_components=6, svd_solver="full")
    vecs = pca.fit_transform(X)

    # Put PCs into a dataframe
    reduced_df = pd.DataFrame(vecs, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"])
    reduced_df["target"] = target.values

    # Flip PC1/PC2 direction to match paper orientation
    reduced_df["PC1"] *= -1
    reduced_df["PC2"] *= -1
    reduced_df.to_csv("CordsforPCA1PCA2")
    # Plot PC1 vs PC2
    plt.figure(figsize=(8, 6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)

    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.show()



def PCAISVD(n_components = 6):

    df = pd.read_csv("Dataset/S3(C)Expression.csv")

    df = df.set_index("ID").T
    df = df.apply(lambda x: np.log2(x+1) if np.issubdtype(x.dtype, np.number) else x)

    
    target = df.index

    imputer = IterativeSVD(rank=8) #Library by fancyimpute
    X_completed = imputer.fit_transform(df.values)

    pca = PCA(n_components=2)
    vecs = pca.fit_transform(X_completed)


    reduced_df = pd.DataFrame(vecs,columns=["PC1", "PC2"])

    reduced_df["target"] = target.values
    reduced_df["PC1"] *= -1
    reduced_df.rename(columns={"PC1": "x", "PC2": "y"}).to_csv("ISVDCords.csv", index=False)



    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.xlabel(f'PC1({pca.explained_variance_ratio_[0]})')
    plt.ylabel(f'PC2({pca.explained_variance_ratio_[1]})')
    plt.show()   



def PCAOALS():
    #Load dataset and exclude ID for loading into array
    df = pd.read_csv("Dataset/S3(C)Expression.csv")
    df = df.apply(lambda x: np.log10(x + 1) if np.issubdtype(x.dtype, np.number) else x)
    
    
    df = df.set_index(df.columns[0]) #tell the dataframe to treat the first column as index and not be calculated into the subsequent calculations
    #Transpose the matrix because PCA usually expects to have rows = samples (IGF Rapa VC etc) and in columns the actual data (similar to S3(A))
    df = df.T
    #After transposing the matrix index is the samples so we use it to later annotate the points in our PCA
    target = df.index
    #We to numpy to treat values that are not specifically called NaN rather _ or .
    dm = df.to_numpy(dtype=float)


    #Calculate the average of every column
    col_means = np.nanmean(dm, axis=0)

    #Find the difference
    dm = dm - col_means

    #Recalculate the values
    n_rows, n_cols = dm.shape
    
    #PCA decomposes the data into Scores(P) and Loadings(T) with a random guess for scores and an empty matrix for loadings
    P = np.random.rand(n_cols,6)
    T = np.zeros((n_rows,6))
    

    #Pold is used later for checking ame for tolerance and max iter
    P_old = np.zeros_like(P)
    tolerance = 1e-6
    max_iter = 100
    

    for iteration in range(max_iter):
            for x in range(n_rows):
                row = dm[x, :]
                mask = ~np.isnan(row) #Map of VALID data (True if number, false if not)
                if not np.any(mask): continue # Skip if row is all NaNs

                P_availiable = P[mask,:] #We create a second list that has only availiable data
                T[x,:] = row[mask] @ np.linalg.pinv(P_availiable.T) #We use least squares (pinv) to find the best loadings for the specific scores

            #Same but for the scores 
            for k in range(n_cols):
                col = dm[:,k]
                mask = ~np.isnan(col)
                if not np.any(mask): continue # Skip if column is all NaNs
                T_availiable = T[mask,:]
                P[k,:] =  np.linalg.pinv(T_availiable) @ col[mask]
        
            #QR decomposition that forces P to be orthogonal
            Q,R = np.linalg.qr(P)
            P = Q
                    

            #Check if the matrix p has stopped changing
            diff = np.linalg.norm(P - P_old)
            if diff < tolerance:
                print(f"Converged at iteration {iteration}")
                break

            P_old = P.copy()

    
    total_var = np.nansum(np.nanvar(dm, axis=0))
    pc_vars = np.var(T,axis=0)
            # Variance explained by each Score column
    explained_vars = pc_vars/total_var
    for i, v in enumerate(explained_vars):
        print(f"PC{i+1} explains {v*100:.2f}% of the variance")

    return P,T,target,explained_vars        


    



def main():
    P,T,target,explVars = PCAOALS()
    reduced_df = pd.DataFrame(T,columns=["PC1","PC2","PC3","PC4","PC5","PC6",])
    reduced_df["target"] = target
    #reduced_df["PC1"] *= -1
    #reduced_df["PC2"] *= -1
    
    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.xlabel(f"PC1 ({explVars[0]*100:.1f}%)")
    plt.ylabel(f"PC1 ({explVars[1]*100:.1f}%)")
    plt.show() 


if __name__ == "__main__":
    main()    