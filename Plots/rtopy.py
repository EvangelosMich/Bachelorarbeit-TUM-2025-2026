import pyreadr as pr
import pandas as pd
import numpy as np
from pyopenms import *
import pymzml as zml
import csv
import os
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind,mannwhitneyu
from urllib.request import urlretrieve
import glob
PPM_TOL = 1.0 #Mass tolerance for example with 100.001 and 100.002 it would be treated as the same molecule
NOISE_THRESHOLD = 5.0 #Any signal below this is treated as electronic static and ignored


DATADIR = "/home/evangelosge84puv/MSV000090865/raw/HILICNEGMZML"
all_mzml_files = glob.glob(f"{DATADIR}/*.mzML")
all_mass_traces = {}
all_peaks = {}
all_featChroms = {}

if not all_mzml_files:
    print("Error no files were found")
else:
    #Initialize mass-trace detection (detect chromatographic traces across scans)
    mtd = MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm",PPM_TOL)# Set m/z tolerance
    mtd_params.setValue("noise_threshold_int",NOISE_THRESHOLD) # Ignore low-intensity noise
    mtd.setParameters(mtd_params)
    #Why: Raw MS data is just a cloud of dots (m/z,rt,intensity). The first step is to connect dots of the same m/z across time
    

    # Initialize elution-peak detection (find chromatographic peaks within traces)
    epd = ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue("width_filtering","fixed")
    epd.setParameters(epd_params)
    #Why:Once mtd has found a line of signals across time (traces),epd looks at that line to see if it foes up and down like a bell curve (what we call a peak) or
    #If its just a flat background line

    # Initialize feature finder (assemble features from detected peaks and traces)
    ffm = FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", "none")
    ffm_params.setValue("remove_single_traces","true")
    ffm_params.setValue("mz_scoring_by_elements","false")
    ffm_params.setValue("report_convex_hulls","true")
    ffm.setParameters(ffm_params)
    #Why:This takes the peaks that were found by epd and bundles them: It calculated the centroid m/z the exact retention time, and the total intensity
    #under the formerly area curve that was found


    fga = FeatureGroupingAlgorithmQT()

    fga_params = fga.getParameters()
    fga_params.setValue("distance_MZ:max_difference", 10.0, "max m/z difference in PPM")
    fga_params.setValue("distance_RT:max_difference", 20.0, "max RT difference in seconds")
    fga.setParameters(fga_params)
    #Setting up the algorithm that will eventually match "Metabolite A" in IGF1.1 with "Metabolite A" in IGF1.2




    
    #Part 2: Processing loop
    #Loading all the raw .mzML files into RAM
    for file in all_mzml_files: 
        #print(file)
        exp = MSExperiment()
        fm_file = FeatureMap()
        MzMLFile().load(file,exp)
        masstraces = []
        masstracessplit = []
        feat_chrom = []
        mtd.run(exp,masstraces,0) #Connecting the dots into traces
        epd.detectPeaks(masstraces,masstracessplit) #Cut traces into specific peaks
        ffm.run(masstracessplit,fm_file,feat_chrom) #Conert peaks into quantified Features
        #Here essentially we create the list of metabolites for each specific file
        all_featChroms[file] = fm_file
        


      


 # TODO extract features like mz,rt,intensity,quality DONE
    fmfileFeatures = []
    for file, fm in all_featChroms.items():   # iterate over dict: (key, value)
        for feature in fm:                    # feature is a Feature object
            fmfileFeatures.append([
                feature.getMZ(),
                feature.getRT(),
                feature.getIntensity(),
                feature.getOverallQuality(),
                feature.getUniqueId(),
                file   # optional: track which sample this came from
            ])
        
    my_columns = [
    'm/z',              # from feature.getMZ()
    'Retention Time',   # from feature.getRT()
    'Intensity',        # from feature.getIntensity()
    'Quality',          # from feature.getOverallQuality()
    'ID',               # from feature.getUniqueId()
    'Filename'          # from file (the key in your dict)
]

    mapList = list(all_featChroms.values())
    consensusMap = ConsensusMap()

    
    for i, (file_path, fmap) in enumerate(all_featChroms.items()):
        fmap.setIdentifier(str(i))

    aligner = MapAlignmentAlgorithmPoseClustering()
    aligner.align(mapList)





    fga.group(mapList, consensusMap)
    

# 1. Create a list to store our rows
consensus_data = []

# 2. Get the filenames to use as column headers
# We use the keys from your all_featChroms dict
filenames = [os.path.basename(f) for f in all_featChroms.keys()]

for feature in consensusMap:
    # Start the row with basic metadata
    row = {
        'mz': feature.getMZ(),
        'rt': feature.getRT(),
        'quality': feature.getQuality()
    }
    
    # Initialize all intensities to 0.0 (or np.nan) for this metabolite
    for name in filenames:
        row[name] = 0.0
        
    # Now fill in the intensities that actually exist for this consensus feature
    # Each 'handle' points to the original file it came from
    for handle in feature.getFeatureList():
        file_index = handle.getMapIndex()
        # Map the intensity to the correct filename column
        row[filenames[file_index]] = handle.getIntensity()
        
    consensus_data.append(row)

# 3. Create the clean DataFrame
df = pd.DataFrame(consensus_data)
df = df.sort_values("mz")
# Save it to check the 1,200 row count
df.to_csv("Consensus_Matrix_Final.csv", index=False)
print(f"Matrix created! Rows: {len(df)}, Columns: {len(df.columns)}")
    





      

 
 # TODO use plot methods to fill in NaN values
  
 
 # TODO perform PCAS
 # TODO Perform T-tests
 # TODO Testing correction with FDR     
    

