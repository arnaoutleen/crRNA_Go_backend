import pandas as pd
import numpy as np
import csv
from Bio.Seq import Seq
from sklearn.model_selection import train_test_split, ShuffleSplit
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, LSTM, MaxPooling1D, Conv1D, Flatten, Bidirectional
import matplotlib.pyplot as plt
from numpy import argmax
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation

def DNA_to_RNA(input_seqn):
	old_seqn=str(input_seqn)
	old_seqn = old_seqn.strip('[')
	old_seqn = old_seqn.strip(']')
	old_seqn = old_seqn.strip("'")
	old_seqn = old_seqn.replace("\r","")
	old_seqn = old_seqn.replace("\n","")
	old_seqn = old_seqn.replace(" ","")
	old_seqn = old_seqn.replace("\\","")
	old_seqn = old_seqn.replace("n","")
	data=old_seqn.replace('T','U').upper()
	return(data)

#function that one-hot encodes the data
def one_hot_this(data):
	alphabet = 'AUGC'
	onehot_encoded_1=[]
	seqn_array_1=[]
	# define a mapping of chars to integers
	char_to_int_1 = dict((c, i) for i, c in enumerate(alphabet))
	int_to_char_1 = dict((i, c) for i, c in enumerate(alphabet))
	# integer encode input data
	integer_encoded_1 = [char_to_int_1[char] for char in data]
	# print(integer_encoded)
	# one hot encode
	onehot_encoded_1 = list()
	for value in integer_encoded_1:
	    letter_1 = [0 for _ in range(len(alphabet))]
	    letter_1[value] = 1
	    onehot_encoded_1.append(letter_1)
	return(onehot_encoded_1)

# # load model
model = load_model('my_model_bidirectional.h5')
model.summary()

## read DNA sequence in 22 bp increment, write to csv, with location and direction
# load sequence
old_seqn=Seq("ACTGACTGCATGCATGCGGGCCAAATTTTGTGAC")
old_rc=old_seqn.reverse_complement()
old_seqn_str=str(old_seqn)
old_rc_str=str(old_rc)

#convert to RNA sequence by replacing T with U
new_seqn_str=old_seqn_str.replace('T','U').upper()
new_rc_str=old_rc_str.replace('T','U').upper()
length_RNA=len(new_rc_str)

steps = 22
#set up loop vars for fwd loop indexing
start_index=0
end_index= steps
temp=[]
temp_2=[]
fwd_strings=[]
fwd_strand=[]
fwd_position=[]
rev_strings=[]
rev_strand=[]
rev_position=[]

#fwd RNA segmentation + writing to array
for i in range(length_RNA):
	if end_index > length_RNA:
		break
	temp=new_seqn_str[start_index:end_index]
	#print(temp)
	fwd_strand.append('+')
	fwd_strings.append(temp)
	fwd_position.append(start_index + 1)
	start_index = start_index + 1
	end_index = end_index + 1

df=pd.DataFrame({'Sequence' : fwd_strings, 'Strand' : fwd_strand,'Position' : fwd_position})
# print(df)

#reset loop vars
start_index=0
end_index= steps
rev_index=length_RNA

#rev RNA segmentation + writing to array
for i in range(length_RNA):
	if end_index > length_RNA:
		break
	temp_2=new_rc_str[start_index:end_index]
	#print(temp_2)
	rev_strings.append(temp_2)
	rev_strand.append('-')
	rev_index = length_RNA - end_index
	rev_position.append(rev_index + 1)
	start_index = start_index + 1
	end_index = end_index + 1

df_2=pd.DataFrame({'Sequence' : rev_strings, 'Strand' : rev_strand,'Position' : rev_position})
new_df=df.append(df_2,ignore_index=True)
# print(new_df)

# # one-hot encode all outputs into a matrix
RNA_inputs=new_df[['Sequence']]

#convert each RNA sequence into array and 
onehot_encoded=[]
seqn_array=[]

RNA_inputs=RNA_inputs.values.tolist()
length=len(RNA_inputs)
# define universe of possible input values
alphabet = 'AUGC'

for i in range(length):
	#convert to RNA sequence by replacing T with U
	old_seqn=str(RNA_inputs[i])
	data = DNA_to_RNA(old_seqn)
	onehot_encoded = one_hot_this(data)
	seqn_array.append(onehot_encoded)

seqn_array = np.asarray(seqn_array)	
seqn_array = np.expand_dims(seqn_array, axis=0)
#input into prediction module
predict_energy=[]
size_of_arg=len(seqn_array)

for k in range(size_of_arg):
	temp_energy = model.predict(seqn_array[k])
	predict_energy.append(temp_energy)
predict_energy = np.asarray(predict_energy)
predict_energy = predict_energy[0, :, :]
# write to dataframe
new_df['Predicted Energy Score'] = predict_energy 
# order by descending energy score
sorted_results = new_df.sort_values('Predicted Energy Score',ascending=False)
sorted_results = sorted_results.reset_index(drop=True)
# print(sorted_results)
# write to csv file
sorted_results.to_csv("results.csv", index=False)
# print(sorted_results.loc[1,"Strand"])

# make a diagram for the top five hits 
gdd = GenomeDiagram.Diagram('Test Diagram')
number = 5
color=[colors.blue, colors.red, colors.green, colors.orange, colors.purple]
# steps
for i in range(number):
	gdt_features_i = gdd.new_track(i+1, greytrack=False)
	gds_features_i = gdt_features_i.new_set()
	temp_start = sorted_results.loc[i+1,"Position"]
	if sorted_results.loc[i+1,"Strand"] == "+":
		temp_end = temp_start + steps
		strand_val = +1
		strand_name = "crRNA{}".format(i+1)
	else:
		temp_end = temp_start + steps
		strand_val = -1
		strand_name = "crRNA{}".format(i+1)
	feature = SeqFeature(FeatureLocation(int(temp_start), int(temp_end), strand=strand_val))
	gds_features_i.add_feature(feature, name=strand_name, label=True, label_size=25, label_angle = 0, sigil="ARROW",color=color[i])	

#draw the diagram and save it
gdd.draw(format='linear', pagesize='A3', fragments=1,
         start=0, end=length_RNA)
gdd.write("GD.jpg", "jpg")