import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load and prepare data (assuming your initial loading is correct)
def plotSpectrum(x,start_time,end_time):
    df = pd.read_csv("Dataset/S3(B)Feature.csv")
    dfImportant = df[['General|All|rtmed','General|All|mzmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']]
    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

    # Use the full, unsliced arrays once
    full_time_array = np.array(dfImportant['General|All|rtmed'])
    full_intensity_array = np.array(dfImportant[x])
    print(dfImportant.head())
    # Initialize the window range
   

    # The upper limit of your desired spectrum is 526, so the loop continues
    # as long as the start_time is below the overall max range.
        
        # 1. Create a FRESH mask based on the full array in EACH iteration
    mask_batch = (full_time_array >= start_time) & (full_time_array <= end_time)
        
        # 2. Slice BOTH arrays using the same mask
    time_batch = full_time_array[mask_batch]
    intensity_batch = full_intensity_array[mask_batch]
    

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
        ax.plot(time_batch,intensity_batch,label=f'{x} intensity')
        ax.scatter(x_maximus,y_maximus,colorizer='r', label = 'local maxima')
        ax.set_xlabel("Retention Time (rtmed)")
        ax.set_title(f'{x} Spectrum({start_time} - {end_time})')
        ax.legend()

    return fig     

