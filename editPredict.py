import h5py 
from keras.models import model_from_json
import random as ran
import numpy as np
import keras
from numpy import array
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--txt", help = "input txt file", required=True)
parser.add_argument("-c", "--json", help = "input txt file", required=True)
parser.add_argument("-w", "--h5", help = "input txt file", required=True)
args = parser.parse_args()

np.set_printoptions(threshold=np.inf)
model = model_from_json(open(args.json).read())
model.load_weights(args.h5)

label_encoder = LabelEncoder()
onehot_encoder = OneHotEncoder(sparse=False)
alphabet = 'ACGT'
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))


with open (args.txt) as tf1:
	
  for line in tf1: 		
  		line =line.upper()
		line=line.strip('\n')
		line= line.split('	')
		if 'N' in line[-1]:
			continue
		
		else:
			l = len(line[-1])
			values=array(list(line[-1]))
   	 	integer_encoded = [char_to_int[char] for char in values]
		onehot_encoded = list()
		for value in integer_encoded:
			letter = [0 for _ in range(len(alphabet))]
			letter[value] = 1
			onehot_encoded.append(letter)
		onehot_encoded = array(onehot_encoded)
		onehot_encoded=onehot_encoded.reshape(1,l,4,1)
		result = model.predict(onehot_encoded)
		result1= model.predict_classes(onehot_encoded)
		print result, result1
		
		


		
		
		


