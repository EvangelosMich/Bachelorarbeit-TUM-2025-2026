import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Load and prepare data (assuming your initial loading is correct)
def plotSpectrum(substance,timeStamp):
    df = pd.read_csv("Dataset/S3(B)Feature.csv")
    dfImportant = df[['General|All|rtmed','General|All|mzmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']] #time mass and mean peak intensities for each indvidual substance excluding quality control
    dfImportant = dfImportant.sort_values(by=['General|All|rtmed']) #sort by time

    # Use the full, unsliced arrays once
    full_time_array = np.array(dfImportant['General|All|rtmed'])
    full_intensity_array = np.array(dfImportant[substance])
    full_mass_array = np.array(dfImportant['General|All|mzmed'])
    #print(dfImportant.head())
    
   

    
    
        
        # 1. Create a mask to keep the times in intervals of 10 seconds
    mask_batch = (full_time_array >= timeStamp - 5) & (full_time_array <= timeStamp + 5)
        
        # 2. Slice BOTH arrays using the same mask
    time_batch = full_time_array[mask_batch] 
    intensity_batch = full_intensity_array[mask_batch]
    mass_batch = full_mass_array[mask_batch]
    

    fig,ax = plt.subplots()
    # Check if the batch is empty before plotting
    if len(time_batch) > 0:
        
        # --- Peak Detection Logic (Modified to use current batch data) ---
        x_maximus = []
        y_maximus = []
        
        # We search for local maxima (peaks) within this small batch
        # Ensure array is long enough for the search window (i-2 to i+2)
        if len(intensity_batch) >= 5:
            for i in range(2, len(intensity_batch) - 2):
                current_i = intensity_batch[i]
                if (current_i > intensity_batch[i-1] and 
                    current_i > intensity_batch[i-2] and
                    current_i > intensity_batch[i+1] and 
                    current_i > intensity_batch[i+2]):
                    
                    y_maximus.append(current_i)
                    x_maximus.append(time_batch[i])

        # --- Plotting ---
        ax.plot(mass_batch,intensity_batch,label=f'{substance} intensity')
        ax.scatter(x_maximus,y_maximus,colorizer='r', label = 'local maxima')
        ax.set_xlabel("Retention Time (rtmed)")
        ax.set_title(f'{substance} Spectrum({timeStamp-5} - {timeStamp+5})')
        ax.legend()

    return fig     


def main():
    fig = (plotSpectrum("Stats|Mean|IGF",255))
    plt.show()
    


if __name__ == "__main__":
    main()    

