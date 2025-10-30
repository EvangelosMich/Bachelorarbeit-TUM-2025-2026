import streamlit as st
import plot 
import pandas as pd
import numpy as np


st.markdown("Peak graphing depending on value")
values = st.slider("Time values",min_value=38,max_value=520)

options = st.multiselect(
    "What substance do you wanna plot?",
    ["IGF-1","Rapamycin","Vehicle control"],
    max_selections=1,
    accept_new_options=False
)
st.write(options)
st.write(values)
if options == ["IGF-1"]:
    st.write("Im inside IGF1")
    fig= plot.plotSpectrum('Stats|Mean|IGF',values-5,values+5)
    st.pyplot(fig)
elif options == ["Rapamycin"]:
    st.write("Im inside Rapa")

    fig= plot.plotSpectrum('Stats|Mean|Rapamycin',values-5,values+5)
    st.pyplot(fig)  
elif options == ["Vehicle control"]:
    st.write("Im inside Control")

    fig= plot.plotSpectrum('Stats|Mean|Control',values-5,values+5)
    st.pyplot(fig)