#!/usr/bin/env python
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
from tensorflow.keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization, ZeroPadding2D, Activation
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D, AveragePooling2D
from tensorflow.keras.layers import concatenate
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1, L2, L1L2
#Set seed
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
        Y_part = np.log(np.genfromtxt(Y_file))
        Y_aggr.append(Y_part)
    Y_aggr = np.array(Y_aggr)
    Y_aggr = np.concatenate(Y_aggr)
    return(Y_aggr)

#CNN regression: inference of branch lengths
def build_CNN_brl(X_train,Y_train,conv_pool_n,filter_n,droput_rates,batch_sizes):
    
    # Length of MSA, i.e. number of sites 
    Aln_length = X_train.shape[2]
    # Number of taxa in MSA
    Ntaxa = X_train.shape[1]
    #Number of branches of a tree
    N_branch = Y_train.shape[1]
    
    
    #1. MSA CNN branch 1 
    #Hyperparameters
    #Hight (horizontal)
    conv_x=[Ntaxa,1,1,1,1,1,1,1]
    #Width (vertical)
    conv_y=[1,2,2,2,2,2,2,2]
    pool=[1,4,4,4,2,2,2,1]
    filter_s=[1024,1024,128,128,128,128,128,128]

    print(f"N convolutional layers: {conv_pool_n}\nN fileters: {filter_n}\nDropout rate: {droput_rates}\nBatch size: {batch_sizes}")
   
    # CNN Arhitecture
    visible_msa = Input(shape=(Ntaxa,Aln_length,1))
    x = visible_msa
    for l in list(range(0,conv_pool_n)):
        x = ZeroPadding2D(padding=((0, 0), (0,conv_y[l]-1)))(x)        
        x = Conv2D(filters=filter_s[l], kernel_size=(conv_x[l], conv_y[l]), strides=1,activation='relu')(x)
        x = Dropout(rate=droput_rates)(x)
        x = BatchNormalization()(x)
        x = AveragePooling2D(pool_size=(1,pool[l]))(x)
        x = Dropout(rate=droput_rates)(x)
    output_msa = Flatten()(x)
    
    
    hidden1 = Dense(1000,activation='relu')(output_msa)
    drop1 = Dropout(rate=droput_rates)(hidden1)
    output = Dense(N_branch, activation='linear')(drop1)
    

    model_cnn = Model(inputs=visible_msa, outputs=output)
    model_cnn.compile(loss='mean_squared_error',optimizer='adam',metrics=['mae','mse'])
    
    #Print model
    print(model_cnn.summary())
   
    #Model stopping criteria
    callback1=EarlyStopping(monitor='val_loss', min_delta=0.001, patience=20, verbose=1, mode='auto')
    callback2=ModelCheckpoint('best_weights_cnn', monitor='val_loss', verbose=0, save_best_only=True, save_weights_only=False, mode='auto', save_freq='epoch')    

    model_cnn.fit(x=X_train,y=Y_train,batch_size=batch_sizes,callbacks=[callback1,callback2],epochs=400,verbose=1,shuffle=True,validation_split=0.1)
    return(model_cnn)
 

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
    print(f"\nConcatenating {args.LABTRAIN} Y datasets in order")
    Y_train = aggregate_Yinput(args.LABTRAIN)
    
    print("\n==========Reading testing data==========")
    print(f"\nConcatenating {args.TEST} X datasets in order")
    X_test = aggregate_Xinput(args.TEST)
    print(f"N taxa: {X_test.shape[1]}\nAlignment length: {X_test.shape[2]}\nN alignments for testing: {X_test.shape[0]}")
    print(f"\nConcatenating {args.LABTEST} Y datasets in order")
    Y_test = aggregate_Yinput(args.LABTEST)


    #Regression BL
    #Run model
    model_cnn_reg=build_CNN_brl(X_train=X_train,Y_train=Y_train,conv_pool_n=6,filter_n=100,droput_rates=0.15,batch_sizes=100)
    
    #Evaluate model
    evals_reg = model_cnn_reg.evaluate(X_test,Y_test,batch_size=100, verbose=1, steps=None)
    bls = model_cnn_reg.predict(X_test,batch_size=100, verbose=1, steps=None)
    np.savetxt("brls.evaluated.cnn.txt",evals_reg,fmt='%f')
    np.savetxt("brls.predicted.cnn.txt",np.exp(bls),fmt='%f')
    
    #Saving model
    print("\nSaving keras trained model")
    model_cnn_reg.save("keras_model_cnn.h5")
    tf.keras.utils.plot_model(model_cnn_reg, to_file='model_cnn.png', show_shapes=True)

    
if __name__ == "__main__":
    main()