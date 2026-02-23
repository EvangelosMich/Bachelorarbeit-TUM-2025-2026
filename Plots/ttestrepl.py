import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def randomTtest():
    #df = pd.read_csv("Dataset/processedFeatureInfo.csv")
    df = pd.read_csv("Dataset/S3(B)Feature.csv")



    plt.figure(figsize=(8,6))
    plt.scatter(-df["ttest|IGF_vs_Control|mean.diff"], df["ttest|IGF_vs_Control|log.fdr"], s=50)
    #plt.scatter(-df["ttest.igf_vs_VC.meandiff"], df["ttest.IGF_vs_VC.logpvalue"], s=50)



    plt.xlabel("-MeanDiff")
    plt.ylabel("Log.fdr")
    plt.xlim(-2.0,2.0)
    plt.title("Volcano Plot Analysis of Differential Features Between Rapamycin and Control Conditions on Original Data")
    plt.show()


def volcanoPythonigfvc():
    df = pd.read_csv("Features_Matrix_Python.csv")
    fig,ax = plt.subplots(figsize=(8,6))

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
    volcanoPythonrapavc()
if __name__ == "__main__":
    main()
