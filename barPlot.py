import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plotSpectrum1(substance, timeStamp,comparison_cols =None,ax=None,color="blue"):
    # Load and sort the relevant columns
    df = pd.read_csv("Dataset/S3(B)Feature.csv")
    dfImportant = df[['General|All|rtmed',"General|All|annot_ms1", 'General|All|mzmed',
                      'Stats|Mean|IGF', 'Stats|Mean|Rapamycin', 'Stats|Mean|Control',"ttest|IGF_vs_Rapamycin|mean.diff","ttest|Control_vs_Rapamycin|mean.diff","ttest|IGF_vs_Control|mean.diff"]]
    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])
    dfFilteredForDiff = df[["General|All|annot_ms1","ttest|IGF_vs_Rapamycin|mean.diff","ttest|Control_vs_Rapamycin|mean.diff","ttest|IGF_vs_Control|mean.diff"]]
    condition1 = (abs(dfFilteredForDiff["ttest|IGF_vs_Rapamycin|mean.diff"]) >= 0.5)
    condition2 = (abs(dfFilteredForDiff["ttest|Control_vs_Rapamycin|mean.diff"]) >= 0.5)
    condition3 = (abs(dfFilteredForDiff["ttest|IGF_vs_Control|mean.diff"])>= 0.5)
    
    if comparison_cols != None and len(comparison_cols)==1:
         dfFilteredForDiff = dfFilteredForDiff[["General|All|annot_ms1",comparison_cols[0]]]
         maskForDiff = (abs(dfFilteredForDiff[comparison_cols[0]])>=0.5)
         dfFilteredForDiff = dfFilteredForDiff[maskForDiff]
    elif comparison_cols != None and len(comparison_cols)==2:
         dfFilteredForDiff = dfFilteredForDiff[["General|All|annot_ms1",comparison_cols[0],comparison_cols[1]]]
         maskForDiff = (abs(dfFilteredForDiff[comparison_cols[0]])>=0.5) |  (abs(dfFilteredForDiff[comparison_cols[1]])>=0.5)
    else:          
        dfFilteredForDiff = dfFilteredForDiff[condition1 | condition2 | condition3]
        maskForDiff = condition1 | condition2 | condition3
    print(dfFilteredForDiff.head())
    print(maskForDiff)    



    pd.set_option('display.max_columns', None)
    # Extract arrays
    full_time_array = np.array(dfImportant['General|All|rtmed'])
    full_intensity_array = np.array(dfImportant[substance])
    full_mass_array = np.array(dfImportant['General|All|mzmed'])

    
    # Select only the data around the given timestamp
    mask_batch = (full_time_array <= timeStamp + 5) & (full_time_array >= timeStamp -5)
    #time_batch = full_time_array[mask_batch]
    mass_batch = full_mass_array[mask_batch]
    intensity_batch = full_intensity_array[mask_batch]
    df_batch_data = dfImportant[mask_batch]

    df_batch_annotate = df_batch_data[maskForDiff]
    df_batch_annotate.to_csv("FilteredCsv.csv")

    # Create pseudo-spectrum plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = plt.gcf() 
    ax.set_xlim(0,1100)
    ax.set_ylim(0,7)
    ax.vlines(mass_batch, 0, intensity_batch, color=color, linewidth=1)

    for index,row in df_batch_annotate.iterrows():
        mz = row["General|All|mzmed"]   
        intensity = row[substance]
        raw_label = row['General|All|annot_ms1']
        label = str(raw_label)
        if ';' in label:
                label = label.split(',')[0].strip()
        
        if pd.notna(label) and label != '':
            ax.text(
            mz, 
            intensity + 0.05, # Place label slightly above the bar
            label, 
            fontsize=5, 
            rotation=0, # Rotate for better readability
            ha='center'  # Center the text above the mz value
        )


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
    plotSpectrum1("Stats|Mean|IGF", 38,comparison_cols=["ttest|IGF_vs_Rapamycin|mean.diff","ttest|IGF_vs_Control|mean.diff"])

if __name__ == "__main__":
    main()
