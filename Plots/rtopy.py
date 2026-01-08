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
PPM_TOL = 1.0
NOISE_THRESHOLD = 5.0


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
    

    # Initialize elution-peak detection (find chromatographic peaks within traces)

    epd = ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue("width_filtering","fixed")
    epd.setParameters(epd_params)


    # Initialize feature finder (assemble features from detected peaks and traces)
    ffm = FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", "none")
    ffm_params.setValue("remove_single_traces","true")
    ffm_params.setValue("mz_scoring_by_elements","false")
    ffm_params.setValue("report_convex_hulls","true")
    ffm.setParameters(ffm_params)
    counter = 0


    fga = FeatureGroupingAlgorithmQT()

    fga_params = fga.getParameters()
    fga_params.setValue("distance_MZ:max_difference", 10.0, "max m/z difference in PPM")
    fga_params.setValue("distance_RT:max_difference", 20.0, "max RT difference in seconds")
    fga.setParameters(fga_params)




    
    
    for file in all_mzml_files: 
        #print(file)
        exp = MSExperiment()
        fm_file = FeatureMap()
        MzMLFile().load(file,exp)
        masstraces = []
        masstracessplit = []
        feat_chrom = []
        mtd.run(exp,masstraces,0)
        epd.detectPeaks(masstraces,masstracessplit)
        ffm.run(masstracessplit,fm_file,feat_chrom)

        all_mass_traces[file] = masstraces
        all_peaks[file] = masstracessplit
        all_featChroms[file] = fm_file
        

        # print(f"Found {len(masstraces)} mass traces")        
        # print(f"Found{len(masstracessplit)} peaks")
        # print(f"Found{fm_file.size()} fms")
    
        # print(f"Just checking something {fm_file[0].getMZ()}")
      


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
    





    #print(fmfileFeatures)        

 
 # TODO use plot methods to fill in NaN values
  
 
 # TODO perform PCAS
 # TODO Perform T-tests
 # TODO Testing correction with FDR     
    

