import pyreadr as pr
import pandas as pd
import numpy as np
from pyopenms import *
import pymzml as zml
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind,mannwhitneyu
from urllib.request import urlretrieve
import glob
PPM_TOL = 15.0
NOISE_THRESHOLD = 50.0


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
        # fm_file.get_df().to_csv(f"out{counter}")
        # counter+=1


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
    print(fmfileFeatures)        
 # TODO merge all featureMaps into one matrix DONEish

 
 # TODO use plot methods to fill in NaN values
  
 
 # TODO perform PCAS
 # TODO Perform T-tests
 # TODO Testing correction with FDR     
    

