import streamlit as st
import matplotlib.pyplot as plt
import barPlot as bp

st.markdown("Peak graphing depending on value")
values = st.slider("Time values",min_value=38,max_value=520) #manually inserting time values can be done with variables

options = st.multiselect( #what Substance wants to be used.
    "What substance do you wanna plot?",
    ["Stats|Mean|IGF","Stats|Mean|Rapamycin","Stats|Mean|Control"],
    max_selections=3,
    accept_new_options=False
)

if len(options) == 1:
    substance = options[0].split('|')[-1] #eg IGF
    if substance == "IGF":
        compare = [f"ttest|{substance}_vs_Rapamycin|mean.diff",f"ttest|{substance}_vs_Control|mean.diff"]
    elif substance == "Rapamycin":
        compare = [f"ttest|IGF_vs_{substance}|mean.diff",f"ttest|Control_vs_{substance}|mean.diff"]
    else:
        compare = [f"ttest{substance}_vs_Rapamycin|mean.diff",f"ttest|IGF_vs_{substance}|mean.diff"]
    


    fig = bp.plotSpectrum1(options[0],values,compare,color='blue')
    st.pyplot(fig)

elif len(options) == 2:
    fig,ax = plt.subplots(figsize = (10,6))
    substance1 = options[0].split('|')[-1]#IGF
    substance2 = options[1].split('|')[-1]#Rapamycin
    compare = f"ttest|{substance1}_vs_{substance2}|mean.diff"
    bp.plotSpectrum1(options[0],values,comparison_cols=[compare],ax=ax,color= "blue")
    
    bp.plotSpectrum1(options[1],values,comparison_cols=[compare],ax=ax,color="red")
    ax.legend(options)
    st.pyplot(fig)

elif len(options)==3:
    fig,ax = plt.subplots(figsize = (7,4))
    bp.plotSpectrum1(options[0],values,ax=ax,color= "blue")
    bp.plotSpectrum1(options[1],values,ax=ax,color="red")
    bp.plotSpectrum1(options[2],values,ax=ax,color='green')
    ax.legend(options)
    st.pyplot(fig)   