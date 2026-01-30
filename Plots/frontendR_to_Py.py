import streamlit as st
import pandas as pd
import os
import rtopy
import vennDiagExpression as ve
import PCA

if "results_ready" not in st.session_state:
    st.session_state.results_ready = False
    st.session_state.df_features = None
    st.session_state.df_expression = None
uploaded_files = st.file_uploader("Choose your files", type=["mzML"],accept_multiple_files=True)


if uploaded_files:
    st.success(f"Loaded {len(uploaded_files)} files")

    if st.button("Run Analysis on Raw metabolomics Data",):
        if uploaded_files:
            resultsFeatures,resultsExpression = rtopy.run_metabolomics_pipeline(uploaded_files)
            st.session_state.df_features = resultsFeatures
            st.session_state.df_expression = resultsExpression
            st.session_state.results_ready = True
        else:
            st.warning("Please upload appropriate Files first")    
            
    if st.session_state.results_ready:
        tab2,tab1 = st.tabs(["Features Matrix","Expression Matrix"])
        
        with tab1:
            st.header("Expressions Table")
            st.dataframe(st.session_state.df_expression)

            if st.button("Normalize according to the original study"):
                
                processedCSV = ve.normalization("Expression_Matrix_Python.csv")
                st.pyplot(PCA.normalPCA(processedCSV))   
        with tab2:
            st.header("Features Table")
            st.dataframe(st.session_state.df_features)

 




