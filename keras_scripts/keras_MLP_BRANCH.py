import tensorflow as tf
import time
import pandas
from itertools import product
import sys, argparse, os
import numpy as np
from math import log, ceil
from scipy.stats import multinomial, chi2, bayes_mvs
from math import factorial
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization, ZeroPadding2D
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D, AveragePooling2D
from tensorflow.keras.layers import concatenate
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.optimizers import Adam

seed = np.random.randint(0,2**32 - 1,1)
np.random.seed(seed)
tf.random.set_seed(seed)
print(f"Using random seed {seed}")

def aggregate_Xinput(X_files_in):
    X_aggr = []
    for X_file in X_files_in:
        X_part = np.load(X_file)
        X_aggr.append(X_part)
    X_aggr = np.array(X_aggr)    
    X_aggr = np.concatenate(X_aggr)
    return(X_aggr)

def aggregate_Yinput(Y_files_in):
    Y_aggr = []
    for Y_file in Y_files_in:
        #Log transforming branch lengths to avoid negative values
        Y_part = np.log(pandas.read_csv(Y_file, delim_whitespace=True, header=None))
        Y_aggr.append(Y_part)
    Y_aggr = pandas.concat(Y_aggr)
    return(Y_aggr)


def get_site_patterns(array_in):
    N_alns = array_in.shape[0]
    N_taxa = array_in.shape[1]
    character_alphabet = ''.join(map(str,np.unique(array_in)))
    site_patterns = list(product(character_alphabet,repeat = N_taxa))
    site_patterns = list(map(''.join,site_patterns))
    dtype = [tuple([i,"f8"]) for i in site_patterns]
    struc_array = np.zeros(N_alns, dtype=dtype)
    #dic_patterns = dict.fromkeys(site_patterns,np.repeat(0,N_alns))
    for i in range(N_alns):
        aln = array_in[i,:,:]
        Aln_length = aln.shape[1] 
        for s in range(Aln_length):     
            array_key = ''.join((map(str,aln[:,s].flatten())))
            struc_array[array_key][i]+=1/Aln_length  
    return(struc_array.view((float,len(struc_array.dtype.names))))     


#Regression: inference of branch lengths of unrooted trees from site patterns  
def build_MLP_brl(X_train,Y_train,droput_rates,batch_sizes):
    
    # Length of feature vector, i.e. number of patterns 
    N_patterns=X_train.shape[1]
    
    #Number of branches of a tree
    N_branch = Y_train.shape[1]
    
    visible_layer = Input(shape=(N_patterns,))
    hidden1 = Dense(10000,activation='relu')(visible_layer)
    drop1 = Dropout(rate=droput_rates)(hidden1)
    hidden2 = Dense(10000,activation='relu')(drop1)
    drop2 = Dropout(rate=droput_rates)(hidden2)
    hidden3 = Dense(10000,activation='relu')(drop2)
    drop3 = Dropout(rate=droput_rates)(hidden3)
    output = Dense(N_branch, activation='linear')(drop3)
    

    model_mlp = Model(inputs=visible_layer, outputs=output)
    model_mlp.compile(loss='mean_squared_error',optimizer='adam',metrics=['mae','mse'])
    
    #Print model
    print(model_mlp.summary())
   
    #Model stopping criteria
    callback1=EarlyStopping(monitor='val_loss', min_delta=0.001, patience=10, verbose=1, mode='auto')
    callback2=ModelCheckpoint('best_weights_mlp', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto', save_freq='epoch')    

    model_mlp.fit(x=X_train,y=Y_train,batch_size=batch_sizes,callbacks=[callback1,callback2],epochs=100,verbose=1,shuffle=True,validation_split=0.1)
    return(model_mlp)


def main():
    parser = argparse.ArgumentParser(description='Keras run')
    parser.add_argument( '--tr', help = "Training MSAs dataset in npy",nargs='+',dest='TRAIN')
    parser.add_argument( '--te', help = "Test MSAs dataset in npy",nargs='+',dest='TEST')
    parser.add_argument( '--trl', help = "Training branch lengths in csv",nargs='+',dest='LABTRAIN')
    parser.add_argument( '--tel', help = "Test branch lengths in csv",nargs='+',dest='LABTEST')
    
    args = parser.parse_args()
    
    #Read inputs
    print("\n==========Reading training data==========")
    print(f"\nConcatenating {args.TRAIN} X datasets in order")
    X_train = aggregate_Xinput(args.TRAIN)
    print(f"N taxa: {X_train.shape[1]}\nAlignment length: {X_train.shape[2]}\nN alignments for training: {X_train.shape[0]}")
    print("\nExtracting site patterns for TRAIN data, it may take a while")
    X_train = get_site_patterns(X_train) 
    Y_train = aggregate_Yinput(args.LABTRAIN)
    
    print("\n==========Reading testing data==========")
    print(f"\nConcatenating {args.TEST} X datasets in order")
    X_test = aggregate_Xinput(args.TEST)
    print(f"N taxa: {X_test.shape[1]}\nAlignment length: {X_test.shape[2]}\nN alignments for testing: {X_test.shape[0]}")
    print("\nExtracting site patterns for TEST data, it may take a while")
    X_test = get_site_patterns(X_test) 
    Y_test = aggregate_Yinput(args.LABTEST)
 

    #Regression BL
    #Run model
    model_mlp_reg=build_MLP_brl(X_train=X_train,Y_train=Y_train,droput_rates=0.15,batch_sizes=100)
    #Evaluate model
    print("Evaluating with best weights")
    evals_reg = model_mlp_reg.evaluate(X_test,Y_test,batch_size=100, verbose=1, steps=None)
    bls = model_mlp_reg.predict(X_test,batch_size=100, verbose=1, steps=None)
    np.savetxt("brls.evaluated.mlp.txt",evals_reg,fmt='%f')
    np.savetxt("brls.predicted.mlp.txt",np.exp(bls),fmt='%f')
    #Saving model
    print("\nSaving keras trained model")
    model_mlp_reg.save("keras_model_mlp.h5")
    tf.keras.utils.plot_model(model_mlp_reg, to_file='model_mlp.png', show_shapes=True)
       
    
if __name__ == "__main__":
    main() 