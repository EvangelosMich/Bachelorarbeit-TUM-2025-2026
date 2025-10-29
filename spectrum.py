import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



df = pd.read_csv("Dataset/S3(B)Feature.csv")

dfImportant = df[['General|All|rtmed','Stats|Mean|IGF','Stats|Mean|Rapamycin','Stats|Mean|Control']]


dfImportant = dfImportant.sort_values(by=['General|All|rtmed'])

dfIGF = dfImportant[['General|All|rtmed','Stats|Mean|IGF']]
dfRapa = dfImportant[['General|All|rtmed','Stats|Mean|Rapamycin']]
dfControl = dfImportant[['General|All|rtmed','Stats|Mean|Control']]

timeArray = np.array(dfImportant['General|All|rtmed'])
intensityArray = np.array(dfImportant['Stats|Mean|Control'])


data1 = 31
data2 = 51
while data2 < 526:
    mask = (timeArray >= data1) & (timeArray <= data2)
    timeBatch = timeArray[mask]
    intensityBatch = intensityArray[mask]
    data1 += 10
    data2 += 10 

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

            plt.figure()
            plt.plot(timeBatch,intensityBatch)
            plt.scatter(x_maximus,y_maximus,color = 'r')
            plt.legend()
            plt.show()       



