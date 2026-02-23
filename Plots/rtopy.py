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
import tempfile

from statsmodels.stats.multitest import multipletests
PPM_TOL = 15.0 #Mass tolerance for example with 100.001 and 100.002 it would be treated as the same molecule
NOISE_THRESHOLD =25.0 #Any signal below this is treated as electronic static and ignored
MIN_SAMPLES = 3
DATADIR = "/home/evangelosge84puv/MSV000090865/raw/HILICNEGMZML"
MIN_FRAC = 0.5
ref = pd.read_csv("/home/evangelosge84puv/Desktop/Bachelorarbeit-TUM-2025-2026/Plots/hmdb_mini_reference.csv")

all_mzml_files = glob.glob(f"{DATADIR}/*.mzML")
def run_metabolomics_pipeline(uploaded_files):
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
        for uploaded_file in uploaded_files: 
            exp = MSExperiment()
            fm_file = FeatureMap()

            # --- START OF FIX ---
            # We create a temporary file on your hard drive 
            with tempfile.NamedTemporaryFile(delete=False, suffix=".mzML") as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp_path = tmp.name # This is the "Real Path" OpenMS needs
            
            try:
                # Use the physical path instead of the buffer
                MzMLFile().load(tmp_path, exp) 
            finally:
                # Always delete the temp file when done to keep your PC clean
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            # --- END OF FIX ---
            masstraces = []
            masstracessplit = []
            feat_chrom = []
            mtd.run(exp,masstraces,0) #Connecting the dots into traces
            epd.detectPeaks(masstraces,masstracessplit) #Cut traces into specific peaks
            ffm.run(masstracessplit,fm_file,feat_chrom) #Conert peaks into quantified Features
            #Here essentially we create the list of metabolites for each specific file
            all_featChroms[uploaded_file.name] = fm_file
            


        


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
    df_filtered = df_filtered[df_filtered["quality"] > 0.7]

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

    expression_df = expression_df.set_index("ID").T
    fillValue = expression_df.median() - np.log10(2)
    expression_df = expression_df.fillna(fillValue)

    expression_df = expression_df.T
    print(expression_df.head())



    


    igf_mean = expression_df[igf_cols].mean(axis=1)
    rapa_mean = expression_df[rapa_cols].mean(axis=1)
    qc_mean = expression_df[qc_cols].mean(axis=1)
    vc_mean = expression_df[vc_cols].mean(axis=1)

    mean_diff_igf_vs_rapa = expression_df[igf_cols].mean(axis=1) - expression_df[rapa_cols].mean(axis=1)
    mean_diff_igf_vs_vc =  expression_df[igf_cols].mean(axis=1) - expression_df[vc_cols].mean(axis=1)
    mean_diff_rapa_vs_vc =  expression_df[rapa_cols].mean(axis=1) - expression_df[vc_cols].mean(axis=1)

    _, p_valuesIGFRapa = ttest_ind(expression_df[igf_cols],expression_df[rapa_cols],axis=1,nan_policy='omit',equal_var=False)
    _, p_valuesIGFVC = ttest_ind(expression_df[igf_cols],expression_df[vc_cols],axis=1,nan_policy='omit',equal_var=False)
    _, p_valuesRapaVC = ttest_ind(expression_df[rapa_cols],expression_df[vc_cols],axis=1,nan_policy='omit',equal_var=False)


    fdr_igfRapa = multipletests(p_valuesIGFRapa, method ='fdr_bh')
    fdr_igfcontrol = multipletests(p_valuesIGFVC, method= "fdr_bh")
    fdr_rapacontrol = multipletests(p_valuesRapaVC, method= "fdr_bh")






    df_features = df_filtered[metadata_cols]

    df_features["annot_ms1"] = "Unknown"
    df_features["annot_ms1"] = df_features["mz"].apply(lambda x: find_hmdb_names(x, mode = 'neg',ppm=10))


    df_features["IGF.Mean"] = igf_mean.values
    df_features["Rapa.Mean"] = rapa_mean.values
    df_features["QC.Mean"] = qc_mean.values
    df_features["VC.Mean"] = vc_mean.values
    df_features["ttest.igf_vs_Rapa.meandiff"] = mean_diff_igf_vs_rapa.values
    df_features["ttest.igf_vs_VC.meandiff"] = mean_diff_igf_vs_vc.values
    df_features["ttest.rapa_vs_VC.meandiff"] = mean_diff_rapa_vs_vc.values
    df_features["ttest.IGF_vs_Rapa.pvalue"] = p_valuesIGFRapa
    df_features["ttest.IGF_vs_Rapa.logpvalue"] = -np.log10(p_valuesIGFRapa)
    df_features["ttest.IGF_vs_VC.pvalue"] = p_valuesIGFVC
    df_features["ttest.IGF_vs_VC.logpvalue"] = -np.log10(p_valuesIGFVC)
    df_features["ttest.VC_vs_Rapa.pvalue"] = p_valuesRapaVC
    df_features["ttest.VC_vs_Rapa.logpvalue"] = -np.log10(p_valuesRapaVC)
    df_features["ttest.IGF_vs_Rapa.fdr"] = fdr_igfRapa[1]
    df_features["ttest.IGF_vs_Rapa.logfdr"] = -np.log10(fdr_igfRapa[1])
    df_features["ttest.Control_vs_Rapa.fdr"] = fdr_rapacontrol[1]
    df_features["ttest.Control_vs_Rapa.logfdr"] = -np.log10(fdr_rapacontrol[1])
    df_features["ttest.IGF_vs_Control.fdr"] = fdr_igfcontrol[1]
    df_features["ttest.IGF_vs_Control.logfdr"] = -np.log10(fdr_igfcontrol[1])



    df_features.to_csv("Features_Matrix_Python.csv",index=False)

    return df_features,expression_df
        
    
    
    # TODO perform PCAS
    # TODO Perform T-tests
    # TODO Testing correction with FDR     
def find_hmdb_names(obs_mz, mode='pos', ppm=10):
    # STEP A: Calculate Neutral Mass
    # Positive [M+H]+: Neutral = m/z - 1.00727
    # Negative [M-H]-: Neutral = m/z + 1.00727
    if mode == 'pos':
        theoretical_neutral = obs_mz - 1.00727
    else:
        theoretical_neutral = obs_mz + 1.00727
    
    # STEP B: Calculate the ppm window
    # tolerance = (mass * ppm) / 1,000,000
    tol = (theoretical_neutral * ppm) / 1e6
    
    # STEP C: Find matches in the reference
    matches = ref[
        (ref['monoisotopic_mass'] >= theoretical_neutral - tol) & 
        (ref['monoisotopic_mass'] <= theoretical_neutral + tol)
    ]
    
    if not matches.empty:
        # Return all names found for this specific m/z
        return "; ".join(matches['name'].unique())
    else:
        return "Unknown"        

def main():
    # P,T,target,explVars = PCAOALS()
    # reduced_df = pd.DataFrame(T,columns=["PC1","PC2"])
    # reduced_df["target"] = target
    # #reduced_df["PC1"] *= -1
    # #reduced_df["PC2"] *= -1
    
    # plt.figure(figsize=(8,6))
    # plt.scatter(reduced_df["PC1"], reduced_df["PC2"], s=50)


    # for i, txt in enumerate(reduced_df["target"]):
    #     plt.annotate(txt, (reduced_df["PC1"][i], reduced_df["PC2"][i]), fontsize=8)
    #plt.title(f"PCA with Imputation technique of OALS for the Dataset S3(C)Expression")
    # plt.xlabel(f"PC1 ({explVars[0]*100:.1f}%)")
    # plt.ylabel(f"PC1 ({explVars[1]*100:.1f}%)")
    # plt.show() 
    run_metabolomics_pipeline(glob.glob(f"{DATADIR}/*.mzML"))

if __name__ == "__main__":
    main()    