import numpy as np
import pandas as pd

#create sets of energetics values: 6 columns/values, each between 0 and 10,000
#70,000 values in each set
y1 = np.random.randint(10000, size=(70000))
y2 = np.random.randint(10000, size=(70000))
y3 = np.random.randint(10000, size=(70000))
y4 = np.random.randint(10000, size=(70000))
y5 = np.random.randint(10000, size=(70000))
y6 = np.random.randint(10000, size=(70000))

#set up vars for DNA seqn generation
length_array=y1.shape[0]
seqn_array=[]

#create random DNA string for each set of energetics values
for i in range(length_array):
    stringy = np.random.randint(4, size=(1000))
    result_string=str(stringy)
    result_string = result_string.replace("0","A")
    result_string = result_string.replace("1","T")
    result_string = result_string.replace("2","G")
    result_string = result_string.replace("3","C")
    old_seqn = result_string.replace(" ","")
    old_seqn = result_string.replace("\r","")
    old_seqn = old_seqn.strip('[')
    old_seqn = old_seqn.strip(']')
    old_seqn = old_seqn.rstrip()
    seqn_array.append(old_seqn)
x = np.asarray(seqn_array)
pd.set_option('display.max_colwidth', None)

#write dataset (seqn + energetics vals) to a DataFrame and then a .csv file
df=pd.DataFrame({'target seqn' : x, 'Access 8nts' : y1,'Access 16nts' : y2,'Energy A.' : y3,'Sequence A.' : y4,'Self Folding' : y5,'Free End' : y6})
df.to_csv("test_data.csv", index=False)