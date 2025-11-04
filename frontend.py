import streamlit as st
import plot 
import pandas as pd
import numpy as np
import barPlot as bp

st.markdown("Peak graphing depending on value")
values = st.slider("Time values",min_value=38,max_value=520) #manually inserting time values can be done with variables

options = st.multiselect( #what Substance wants to be used.
    "What substance do you wanna plot?",
    ["IGF-1","Rapamycin","Vehicle control"],
    max_selections=1,
    accept_new_options=False
)

if options == ["IGF-1"]:
    st.write("Im inside IGF1")
    fig= bp.plotSpectrum1('Stats|Mean|IGF',values)
    st.pyplot(fig)
elif options == ["Rapamycin"]:
    st.write("Im inside Rapa")

    fig= bp.plotSpectrum1('Stats|Mean|Rapamycin',values)
    st.pyplot(fig)  
elif options == ["Vehicle control"]:
    st.write("Im inside Control")

    fig= bp.plotSpectrum1('Stats|Mean|Control',values)
    st.pyplot(fig)