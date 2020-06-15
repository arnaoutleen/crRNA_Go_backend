import pandas as pd
import numpy as np
import csv
from sklearn.model_selection import train_test_split, ShuffleSplit
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, LSTM, MaxPooling1D, Conv1D, Flatten
import matplotlib.pyplot as plt
from keras.models import load_model
from numpy import argmax

# # Random training data
# x_train = np.random.rand(4026, 5)
# # Random training labels
# y_train = np.random.randint(0, 2, 4026)
# print(x_train.shape, y_train.shape)

#function that converts DNA input to RNA seqn + all uppercase
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


stringy_1 = np.random.randint(4, size=(1000))
result_string_1=str(stringy_1)
result_string_1 = result_string_1.replace("0","A")
result_string_1 = result_string_1.replace("1","T")
result_string_1 = result_string_1.replace("2","G")
result_string_1 = result_string_1.replace("3","C")
data = DNA_to_RNA(result_string_1)

onehot_encoded_1 = one_hot_this(data)
onehot_encoded_1 = np.asarray(onehot_encoded_1)
onehot_encoded_1 = np.expand_dims(onehot_encoded_1, axis=0)