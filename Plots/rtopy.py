import pyreadr as pr
import pandas as pd
import numpy as np
from pyopenms import *
import pymzml as zml
import csv
import os
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind
from urllib.request import urlretrieve
import glob
import PCA as pca
PPM_TOL = 15.0 #Mass tolerance for example with 100.001 and 100.002 it would be treated as the same molecule
NOISE_THRESHOLD =25.0 #Any signal below this is treated as electronic static and ignored
MIN_SAMPLES = 2
DATADIR = "/home/evangelosge84puv/MSV000090865/raw/HILICNEGMZML"
MIN_FRAC = 0.5

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







    fga.group(mapList, consensusMap)

    aligner = MapAlignmentAlgorithmPoseClustering()
    align_params = aligner.getDefaults()
    aligner.setParameters(align_params)

    transformations = []

    for fmap in mapList:
        trafo = TransformationDescription()
        if not fmap.isMetaEmpty():
            aligner.align(fmap, trafo)
            MapAlignmentTransformer.transformRetentionTimes(fmap,trafo)    
    

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
#df.to_csv("Consensus_Matrix_Final.csv", index=False)
#print(f"Matrix created! Rows: {len(df)}, Columns: {len(df.columns)}")

intensity_cols = filenames
metadata_cols = ['mz','rt','quality']
df = df.dropna(thresh=MIN_SAMPLES, subset=intensity_cols)
df = df[df["quality"] > 0.5]
df_metadata = pd.read_csv("Dataset/metadata.csv")
metadata_dic = dict(zip(df_metadata["Filename"],df_metadata["Group"]))


presence = df[intensity_cols].notna()
mask = presence.sum(axis=1) >= MIN_SAMPLES
df_filtered = df[mask]

mask = presence.mean(axis=1) >= MIN_FRAC
df_filtered = df[mask]

df_intensities = df_filtered[intensity_cols].replace(0.0,np.nan)
df_intensities.rename(columns=metadata_dic,inplace=True)
sorted_cols = sorted(df_intensities.columns, key = lambda x: metadata_dic.get(x,x))
df_grouped = df_intensities[sorted_cols]
df_log = np.log10(df_grouped)

expression_df = pd.DataFrame(np.log10(df_grouped), columns=df_grouped.columns)
global_median = expression_df.stack().median()

for col in expression_df.columns:
    medianVal = expression_df[col].median()
    expression_df[col] = expression_df[col] - medianVal + global_median


expression_df.insert(0,'ID',[f"FT{i+1:04d}" for i in range(len(df))])
expression_df.to_csv("Expression_Matrix_Python.csv",index=False)


igf_cols = expression_df.filter(like="IGF1").columns
rapa_cols = expression_df.filter(like="Rapamycin").columns
qc_cols = expression_df.filter(like="QC").columns
vc_cols = expression_df.filter(like="VControl").columns

expression_df[igf_cols] = expression_df.fillna(expression_df[igf_cols].min(axis=1)-np.log10(2))

igf_mean = expression_df[igf_cols].mean(axis=1)
rapa_mean = expression_df[rapa_cols].mean(axis=1)
qc_mean = expression_df[qc_cols].mean(axis=1)
vc_mean = expression_df[vc_cols].mean(axis=1)

print(expression_df.head())

mean_diff_igf_vs_rapa = expression_df[igf_cols].mean(axis=1) - expression_df[rapa_cols].mean(axis=1)
stat, p_val = ttest_ind(expression_df[igf_cols], expression_df[rapa_cols], nan_policy='omit')
print(stat)








df_features = df_filtered[metadata_cols]
df_features.loc[:, [
    "IGF.Mean",
    "Rapa.Mean",
    "QC.Mean",
    "VC.Mean",
    "ttest.igf_vs_rapa.meandiff"
]] = [
    igf_mean,
    rapa_mean,
    qc_mean,
    vc_mean,
    mean_diff_igf_vs_rapa
]



df_features.to_csv("Features_Matrix_Python.csv",index=False)

    
  
 
 # TODO perform PCAS
 # TODO Perform T-tests
 # TODO Testing correction with FDR     
    

