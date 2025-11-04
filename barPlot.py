import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plotSpectrum1(substance, timeStamp):
    # Load and sort the relevant columns
    df = pd.read_csv("Dataset/S3(B)Feature.csv")
    dfImportant = df[['General|All|rtmed', 'General|All|mzmed',
                      'Stats|Mean|IGF', 'Stats|Mean|Rapamycin', 'Stats|Mean|Control']]
    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

    # Extract arrays
    full_time_array = np.array(dfImportant['General|All|rtmed'])
    full_intensity_array = np.array(dfImportant[substance])
    full_mass_array = np.array(dfImportant['General|All|mzmed'])

    dfImportant.to_csv("SortedCSV.csv",index=False)
    
    # Select only the data around the given timestamp
    mask_batch = (full_time_array <= timeStamp + 5) & (full_time_array >= timeStamp -5)
    #time_batch = full_time_array[mask_batch]
    mass_batch = full_mass_array[mask_batch]
    intensity_batch = full_intensity_array[mask_batch]

    # Create pseudo-spectrum plot
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.vlines(mass_batch, 0, intensity_batch, color='steelblue', linewidth=1)
    ax.set_xlabel("m/z")
    ax.set_ylabel(f"Intensity (mean) {substance}")
    ax.set_title(f"Pseudo-spectrum at RT ≈ {timeStamp} s")

    # Optional: make it visually cleaner
    # ax.grid(False)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()
    return fig

def main():
    plotSpectrum1("Stats|Mean|IGF", 252)

if __name__ == "__main__":
    main()
