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
    df = pd.read_csv("Dataset/S3(B)Feature.csv")
    st.write("Im in the IGF-1 Group")
    dfImportant = df[['General|All|rtmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']]


    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

    dfIGF = dfImportant[['General|All|rtmed','Stats|Mean|IGF']]
    #dfRapa = dfImportant[['General|All|rtmed','Stats|Mean|Rapamycin']]
    #dfControl = dfImportant[['General|All|rtmed','Stats|Mean|Control']]

    timeArray = np.array(dfImportant['General|All|rtmed'])
    intensityArray = np.array(dfImportant['Stats|Mean|IGF'])


    data1 = values - 5
    data2 = values + 5

    mask = (timeArray >= data1) & (timeArray <= data2)
    timeBatch = timeArray[mask]
    intensityBatch = intensityArray[mask]

    if len(timeBatch) > 0:
                x_maximus = []
                y_maximus = []
                if len(intensityBatch) >=5:
                    for i in range(2,len(intensityBatch)-2):
                            curr_i = intensityBatch[i]
                            if (curr_i > intensityBatch[i-1] and
                                curr_i > intensityBatch[i-2] and
                                curr_i > intensityBatch[i+1] and
                                curr_i > intensityBatch[i+1]):
                                y_maximus.append(curr_i)
                                x_maximus.append(timeBatch[i])

                    plot.plt.figure()
                    plot.plt.plot(timeBatch,intensityBatch)
                    plot.plt.scatter(x_maximus,y_maximus,color = 'r')
                    plot.plt.legend()
                    st.pyplot(plot.plt.show())       
if options == ["Rapamycin"]:
    st.write("Im in the Rapamycin group")
    df = pd.read_csv("Dataset/S3(B)Feature.csv")

    dfImportant = df[['General|All|rtmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']]


    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

    dfIGF = dfImportant[['General|All|rtmed','Stats|Mean|IGF']]
    #dfRapa = dfImportant[['General|All|rtmed','Stats|Mean|Rapamycin']]
    #dfControl = dfImportant[['General|All|rtmed','Stats|Mean|Control']]

    timeArray = np.array(dfImportant['General|All|rtmed'])
    intensityArray = np.array(dfImportant['Stats|Mean|Rapamycin'])


    data1 = values - 10
    data2 = values + 10

    mask = (timeArray >= data1) & (timeArray <= data2)
    timeBatch = timeArray[mask]
    intensityBatch = intensityArray[mask]

    if len(timeBatch) > 0:
                x_maximus = []
                y_maximus = []
                if len(intensityBatch) >=5:
                    for i in range(2,len(intensityBatch)-2):
                            curr_i = intensityBatch[i]
                            if (curr_i > intensityBatch[i-1] and
                                curr_i > intensityBatch[i-2] and
                                curr_i > intensityBatch[i+1] and
                                curr_i > intensityBatch[i+1]):
                                y_maximus.append(curr_i)
                                x_maximus.append(timeBatch[i])

                    plot.plt.figure()
                    plot.plt.plot(timeBatch,intensityBatch)
                    plot.plt.scatter(x_maximus,y_maximus,color = 'r')
                    plot.plt.legend()
                    st.pyplot(plot.plt.show())
if options == ["Vehicle control"]:
    st.write("IM in the control group")
    df = pd.read_csv("Dataset/S3(B)Feature.csv")

    dfImportant = df[['General|All|rtmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']]


    dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

    dfIGF = dfImportant[['General|All|rtmed','Stats|Mean|IGF']]
    #dfRapa = dfImportant[['General|All|rtmed','Stats|Mean|Rapamycin']]
    #dfControl = dfImportant[['General|All|rtmed','Stats|Mean|Control']]

    timeArray = np.array(dfImportant['General|All|rtmed'])
    intensityArray = np.array(dfImportant['Stats|Mean|Control'])


    data1 = values - 10
    data2 = values + 10

    mask = (timeArray >= data1) & (timeArray <= data2)
    timeBatch = timeArray[mask]
    intensityBatch = intensityArray[mask]

    if len(timeBatch) > 0:
                x_maximus = []
                y_maximus = []
                if len(intensityBatch) >=5:
                    for i in range(2,len(intensityBatch)-2):
                            curr_i = intensityBatch[i]
                            if (curr_i > intensityBatch[i-1] and
                                curr_i > intensityBatch[i-2] and
                                curr_i > intensityBatch[i+1] and
                                curr_i > intensityBatch[i+1]):
                                y_maximus.append(curr_i)
                                x_maximus.append(timeBatch[i])

                    plot.plt.figure()
                    plot.plt.plot(timeBatch,intensityBatch)
                    plot.plt.scatter(x_maximus,y_maximus,color = 'r')
                    plot.plt.legend()
                    st.pyplot(plot.plt.show())                                                      

