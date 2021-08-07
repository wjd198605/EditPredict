import sys
import re
import numpy as np
import keras
import random as ran
import h5py
from keras.models import Sequential
from keras.layers import Dense, Dropout,Activation,Flatten
from keras import optimizers
from numpy import array
from numpy import argmax
from keras.utils import to_categorical
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from keras.layers import Conv2D, MaxPooling2D
from keras.layers.normalization import BatchNormalization
from keras import regularizers
from keras.utils import plot_model
from keras import initializers
from sklearn.model_selection import train_test_split
from keras import backend as K
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix
import itertools
import argparse
from argparse import RawTextHelpFormatter
import os

os.environ['KMP_DUPLICATE_LIB_OK']='True'
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--txt", help = "input sequence txt file", required=True)
args = parser.parse_args()
 

label_encoder = LabelEncoder()
onehot_encoder = OneHotEncoder(sparse=False)
dataseq=[]
datala=[]

alphabet = 'ACGT'
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))
with open (args.txt) as f1:
  for line in f1:
    line=line.upper()
    line=line.strip()
    length= len(line[0:-1])
    if "N" not in line and 'M' not in line and 'K' not in line and 'W' not in line and 'R' not in line and 'Y' not in line:
      values=array(list(line[0:-1]))
      integer_encoded = [char_to_int[char] for char in values]
      onehot_encoded = list()
      for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1
        onehot_encoded.append(letter)
      onehot_encoded = array(onehot_encoded)
      dataseq.append(onehot_encoded)
      datala.append(line[-1])
dataseq = np.array(dataseq)
datala = map(int,datala)

img_rows, img_cols =  length,4

x_train, x_test, y_train, y_test = train_test_split(dataseq, datala, test_size=0.2)
print x_train[0].shape
x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols,1)
x_test =x_test.reshape(x_test.shape[0], img_rows, img_cols,1)


batch_size = 16
num_classes = 2
epochs = 20
input_shape = (img_rows, img_cols,1)


print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')


y_train =  keras.utils.to_categorical(y_train, num_classes)
y_test1  =  keras.utils.to_categorical(y_test, num_classes)
model = Sequential()

model.add(Conv2D(128,(8,4),input_shape=input_shape))
model.add(Activation('relu'))
model.add(BatchNormalization())
model.add(Conv2D(128,(8,1)))
model.add(Activation('relu'))
model.add(BatchNormalization())
model.add(MaxPooling2D(pool_size=(2,1)))

model.add(Conv2D(64,(3,1)))
model.add(Activation('relu'))
model.add(BatchNormalization())
model.add(Conv2D(64,(3,1)))
model.add(Activation('relu'))
model.add(BatchNormalization())
model.add(MaxPooling2D(pool_size=(2,1)))
model.add(Dropout(0.5))
 
model.add(Flatten())
model.add(Dense(128,  activation='relu' ))
model.add(Dense(64,  activation='relu'))
model.add(Dense(num_classes, activation='softmax'))
 
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.SGD(lr=0.01, momentum=0.0, decay=0.0, nesterov=False),
              metrics=['accuracy'])
 
model.fit(x_train, y_train,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test1))
score = model.evaluate(x_test, y_test1, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])

json_string =model.to_json()
open ('EditPredict.json','w').write(json_string)
model.save_weights('EditPredict.h5')




