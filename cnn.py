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

#read .csv file that contains all the data
data=pd.read_csv("test_data.csv", sep=',')
#split .csv file into scores and DNA sequences
energetics_scores=data[['Access 8nts','Access 16nts','Energy A.','Sequence A.','Self Folding','Free End']]
sequences=data[['target seqn']]

#sum up energetics scores in each row to one score
for row in energetics_scores:
	rowtotal=0
	for column in row[1:]:
		rowtotal = (energetics_scores.sum(axis=1))
		rowtotal = rowtotal / 60000
		# print(rowtotal)
#print the sum to one numpy array for later use
refined_with_sum=energetics_scores.sum(axis=1).to_numpy()
refined_with_sum=refined_with_sum / 60000

# print(refined_with_sum)
# print(sequences)


#convert each DNA sequence into array and one-hot encode it
seqn_array=[]
sequences=sequences.values.tolist()
length=len(sequences)
# define universe of possible input values
for i in range(length):
	#convert to RNA sequence by replacing T with U
	data = DNA_to_RNA(sequences[i])
	one_hot_encoded = one_hot_this(data)
	seqn_array.append(one_hot_encoded)
print(seqn_array)

# #################################################################################################
# now, insert ML algo


# #split into train and test groups
# #convert np array to tensor
x = np.asarray(seqn_array)
y = np.asarray(refined_with_sum)
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25)
x_train=tf.convert_to_tensor(x_train,dtype=tf.int32)
y_train=tf.convert_to_tensor(y_train,dtype=tf.float32)
x_test=tf.convert_to_tensor(x_test,dtype=tf.int32)
y_test=tf.convert_to_tensor(y_test,dtype=tf.float32)
# print(x_train.shape)

# define model - two 1D Convolution layers, with 20% dropout, flatten then dense
model = Sequential()
model.add(Conv1D(filters=20, kernel_size=9, strides=1, input_shape=(1000,4)))
model.add(Dropout(0.2))
model.add(Conv1D(filters=4, kernel_size=9, input_shape=(1000,4)))
model.add(Dropout(0.2))
model.add(MaxPooling1D(pool_size=2,strides=2))
model.add(Flatten())
model.add(Dense(1, activation='relu'))
#compile model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.summary()
#fit model
model.fit(x_train, y_train, epochs=1, validation_data=(x_test, y_test))
#save model
model.save('my_model_cnn.h5')

# #######################################################################################################

# # save model output loss vals as images

# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
# plt.show()
plt.savefig('loss.png')


#############################################################################################
#make test DNA seqn to verify the model can predict

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
# onehot_encoded_1.reshape((1,1000,4))

#####################################################################\
#try it out
a=model.predict(onehot_encoded_1)
print(a)