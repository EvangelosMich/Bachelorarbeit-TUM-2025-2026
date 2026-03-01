import pandas as pd
from sklearn.decomposition import PCA   
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from fancyimpute import IterativeSVD
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.metrics import mean_squared_error
import csv
import seaborn as sns
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"], # Matches standard LaTeX font
    "text.usetex": False,                     # Set to True if you have TeX installed on your PC
    "axes.labelsize": 12,                     # Size of X and Y labels
    "font.size": 11,                          # Matches your 11pt document size
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (8, 5),                 # Golden ratio-ish for 1\linewidth
    "savefig.dpi": 300,                       # High resolution
    "savefig.bbox": 'tight'                   # Removes unnecessary white margins
})
def normalPCA(csvFile = None):
    # Load data and transpose so rows=samples, cols=features
    if csvFile is None:
        df = pd.read_csv("Dataset/ExpressionProcessed.csv").set_index("ID").T
    else:
        df = pd.read_csv(csvFile).set_index("ID").T
#     for col in df.select_dtypes(include=[object]):
#         df[col] = df[col].replace(r'^\s*$', np.nan, regex=True)
#     df.index = [name.split('_')[0] for name in df.index]
#     nan_percentage = (df.isnull().sum().sum() / df.size) * 100
#     print(nan_percentage)
# # Or, if you want to keep the replicate number (e.g., IGF1_1)
#     df.index = [name.replace('_Met_PB_06.12.20', '') for name in df.index]
#     plt.figure(figsize=(15, 10), dpi=130)

#     ax = plt.axes()
#     my_cmap = sns.color_palette(["#000000", "#FF0000"])
#     sns.heatmap(df.isna().transpose(),cmap=my_cmap, cbar=False, ax=ax, 
#                 xticklabels=df.index, yticklabels=False)    
    
#     plt.title("Missing Values", fontsize=20)
#     plt.xlabel("Experimental Samples", fontsize = 20) 
#     plt.ylabel("Metabolite Features(ID)", fontsize = 20)
    
#     plt.tight_layout()
#     plt.show()
    # Log10 transform (keep values from exploding, avoid log(0) with +1)
    #df = df.apply(lambda x: np.log10(x + 1) if np.issubdtype(x.dtype, np.number) else x)
    df_pre_fill = df.copy()
    # Sample labels for plotting
    target = df.index
    fillValues = df.median() - np.log10(2)

    # Fill missing values + remove zero-variance features
    df = df.fillna(fillValues)
    df = df.loc[:, df.var(ddof=1) > 0]

    # Run PCA
    X = df.to_numpy(dtype=float)
    pca = PCA(n_components=6, svd_solver="full")
    

    vecs = pca.fit_transform(X)
    
    X_reconstructed = pca.inverse_transform(vecs)

    # 2. Identify where we had real data (not the imputed values)
    # We need to know which values in the original 'df' (pre-filling) were not NaN
    original_data_mask = ~df_pre_fill.isna().to_numpy() 

    # 3. Calculate the error only for those 'real' points
    # We can use (Original - Reconstructed)^2
    error = (X[original_data_mask] - X_reconstructed[original_data_mask])**2
    mse = np.mean(error)

    # Put PCs into a dataframe
    reduced_df = pd.DataFrame(vecs, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"])
    reduced_df["target"] = target.values

    # Flip PC1/PC2 direction to match paper orientation
    reduced_df["PC1"] *= -1
    #reduced_df["PC2"] *= -1
    #reduced_df.to_csv("CordsforPCA1PCA2")
    # Plot PC1 vs PC2
    fig,ax = plt.subplots(figsize=(15, 10), dpi=130)
    scat = ax.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)
    
    for i, txt in enumerate(reduced_df["target"]):
        ax.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=11)

    points = scat.get_offsets()
    x_data = np.array(points[:,0])
    y_data = np.array(points[:,1])
    print(f"Reconstruction MSE: {mse}")
    ax.set_title(f"PCA with Imputation technique of median - log_10(2) for the Python Produced Dataset")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.show()
    return fig



def PCAISVD(n_components = 6):

    df = pd.read_csv("Dataset/S3(C)Expression.csv")

    df = df.set_index("ID").T
    df = df.apply(lambda x: np.log2(x+1) if np.issubdtype(x.dtype, np.number) else x)

    
    target = df.index

    imputer = IterativeSVD(rank=8) #Library by fancyimpute
    X_completed = imputer.fit_transform(df.values)

    pca = PCA(n_components=6)
    vecs = pca.fit_transform(X_completed)
    X_reconstructed = pca.inverse_transform(vecs)
    original_data_mask = df.isna().to_numpy() 

    # 3. Calculate the error only for those 'real' points
    # We can use (Original - Reconstructed)^2
    error = (X_completed[original_data_mask] - X_reconstructed[original_data_mask])**2
    mse = np.mean(error)
    print(f"Reconstruction MSE: {mse}")

    reduced_df = pd.DataFrame(vecs,columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"])

    reduced_df["target"] = target.values
    reduced_df["PC1"] *= -1
    reduced_df.rename(columns={"PC1": "x", "PC2": "y"}).to_csv("ISVDCords.csv", index=False)



    plt.figure(figsize=(8,6))
    plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


    for i, txt in enumerate(reduced_df["target"]):
        plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)

    plt.title(f"PCA with Imputation technique of ISVD for the Dataset S3(C)Expression")
    plt.xlabel(f'PC1({pca.explained_variance_ratio_[0]})')
    plt.ylabel(f'PC2({pca.explained_variance_ratio_[1]})')
    plt.show()   



def PCA_OALS(dm_centered, n_pcs=2, max_iter=100, tol=1e-6, lmbda=0.01):
    n_rows, n_cols = dm_centered.shape
    
    dm_filled = np.nan_to_num(dm_centered, nan=0.0)
    U, s, Vt = np.linalg.svd(dm_filled, full_matrices=False)
    P = Vt[:n_pcs, :].T 
    T = U[:, :n_pcs] * s[:n_pcs]

    for iteration in range(max_iter):
        P_old = P.copy()
        
        for i in range(n_rows):
            mask = ~np.isnan(dm_centered[i, :])
            P_avail = P[mask, :]
            rhs = P_avail.T @ dm_centered[i, mask]
            lhs = P_avail.T @ P_avail + lmbda * np.eye(n_pcs)
            T[i, :] = np.linalg.solve(lhs, rhs)

        for j in range(n_cols):
            mask = ~np.isnan(dm_centered[:, j])
            T_avail = T[mask, :]
            rhs = T_avail.T @ dm_centered[mask, j]
            lhs = T_avail.T @ T_avail + lmbda * np.eye(n_pcs)
            P[j, :] = np.linalg.solve(lhs, rhs)

        Q, R = np.linalg.qr(P)
        P = Q
        
        if np.linalg.norm(P - P_old) < tol:
            break

    total_ss = np.nansum(dm_centered**2)
    explVars = []
    for i in range(n_pcs):

        pc_variance = np.sum(T[:, i]**2) / total_ss
        explVars.append(pc_variance)

    return P, T, explVars


    
def calculate_mse(dm_centered, T, P):
    # 1. Reconstruct the full matrix from the PCA model
    # Note: If P was returned as (n_pcs, n_genes), use P.T
    dm_recon = T @ P.T 
    
    # 2. Create a mask of available (non-missing) data
    mask = ~np.isnan(dm_centered)
    
    # 3. Calculate error only for observed entries
    errors = dm_centered[mask] - dm_recon[mask]
    
    # 4. Calculate Mean Squared Error
    mse = np.mean(errors**2)
    
    # Optional: Calculate RMSE (often easier to interpret)
    rmse = np.sqrt(mse)
    
    return mse, rmse


def main():
# 1. Load Data
    # df = pd.read_csv("Dataset/S3(C)Expression.csv")
    
    # # 2. Pre-process: Log transform and Transpose
    # # So that rows = samples (targets) and columns = genes (variables)
    # num_cols = df.select_dtypes(include=[np.number]).columns
    # df[num_cols] = np.log10(df[num_cols] + 1)
    # df_transposed = df.set_index(df.columns[0]).T
    
    # # Save the index for your plot labels
    # target = df_transposed.index
    
    # # 3. Create the matrix and CENTER it
    # # We MUST use np.nanmean so the mean isn't 'NaN'
    # dm = df_transposed.to_numpy(dtype=float)
    # col_means = np.nanmean(dm, axis=0)
    # dm_centered = dm - col_means # This is what you pass to the function

    # # 4. CALL THE FUNCTION
    # # Pass dm_centered and specify 2 components
    # P, T, explVars = PCA_OALS(dm_centered, n_pcs=2)

    # # 5. Plotting
    # reduced_df = pd.DataFrame(T, columns=["PC1", "PC2"])
    # reduced_df["target"] = target

    # # Flip PC1 to match standard orientation if needed
    # reduced_df["PC1"] *= -1
    # reduced_df["PC2"] *= -1

    # reduced_df.to_csv("PCA_Results_OALS.csv", index=False)
    # print("PCA coordinates saved to 'PCA_Results_OALS.csv'")
    # plt.figure(figsize=(8,6))
    # plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)
    # print(calculate_mse(dm_centered,T,P))
    # for i, txt in enumerate(reduced_df["target"]):
    #     plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)
    
    # plt.title(f"PCA (O-ALS) for Dataset S3(C)Expression")
    # # explVars are fractions (e.g. 0.45), so multiply by 100
    # plt.xlabel(f"PC1 ({explVars[0]*100:.2f}%)")
    # plt.ylabel(f"PC2 ({explVars[1]*100:.2f}%)")
    # plt.axhline(0, color='grey', lw=1, alpha=0.5)
    # plt.axvline(0, color='grey', lw=1, alpha=0.5)
    # plt.show()
    normalPCA("Dataset/ExpressionProcessed.csv")


if __name__ == "__main__":
    main()    