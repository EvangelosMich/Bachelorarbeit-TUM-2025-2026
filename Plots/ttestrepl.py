import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
def randomTtest():
# Load data
    df = pd.read_csv("Dataset/Features_Aligned_Final.csv")

    fig, ax = plt.subplots()
    
    # Calculate colors based on significance (Standard Volcano Logic)
    # Assuming log.fdr > 1.3 (~p < 0.05) is significant
    colors = ['red' if (abs(x) > 0.5 and y > 1.3) else 'gray' 
              for x, y in zip(df["ttest.Rapamycin_vs_Control.mean.diff"], df["ttest.Rapamycin_vs_Control.log.fdr"])]

    ax.scatter(-df["ttest.Rapamycin_vs_Control.mean.diff"], 
               df["ttest.Rapamycin_vs_Control.log.fdr"], 
               c=colors, alpha=0.6, s=30, edgecolors='none')

    ax.set_xlabel(r"$- \Delta \text{Mean (log}_{10}\text{)}$") # Using LaTeX math notation
    ax.set_ylabel(r"$-\log_{10}(\text{FDR})$")
    ax.set_xlim(-2.0, 2.0)
    ax.axhline(y=1.3, color='black', linestyle='--', linewidth=0.8) # Significance threshold line
    
    ax.set_title("Volcano Plot: Rapamycin vs Control (Original Data)")
    
    # Save as PDF for the thesis
    plt.savefig("RapaControl_Original.pdf")
    plt.show()


def volcanoPythonigfvc():
    df = pd.read_csv("Features_Matrix_Python.csv")
    fig,ax = plt.subplots(figsize=(15,10),dpi = 130)

    scat = ax.scatter(df["ttest.igf_vs_VC.meandiff"], df["ttest.IGF_vs_VC.logpvalue"], s=50)
    ax.set_xlabel("-MeanDiff")
    ax.set_ylabel("Log.pvalue")
    ax.set_xlim(-2.0,2.0)
    ax.set_title("Volcano Plot Analysis of Differential Features Between Rapamycin and Control Conditions on Python Data")
    plt.show()
    return fig
def volcanoPythonrapavc():
    df = pd.read_csv("Features_Matrix_Python.csv")
    fig,ax = plt.subplots(figsize=(8,6))

    scat = ax.scatter(df["ttest.rapa_vs_VC.meandiff"], df["ttest.VC_vs_Rapa.logpvalue"], s=50)
    ax.set_xlabel("-MeanDiff")
    ax.set_ylabel("Log.pvalue")
    ax.set_xlim(-2.0,2.0)
    ax.set_title("Volcano Plot Analysis of Differential Features Between Rapamycin and Control Conditions on Python Data")
    plt.show()
    return fig


def main():
    randomTtest()
if __name__ == "__main__":
    main()
