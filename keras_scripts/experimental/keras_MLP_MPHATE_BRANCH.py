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
from tensorflow import keras
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, Dense, Flatten, Dropout, BatchNormalization, ZeroPadding2D
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D, AveragePooling2D
from tensorflow.keras.layers import concatenate
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.optimizers import Adam
import m_phate
import m_phate.train
import scprep

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


def get_site_patterns(array_in):
    N_alns = array_in.shape[0]
    N_taxa = array_in.shape[1]
    character_alphabet = list(map(str,np.unique(array_in)))
    character_alphabet = [item for item in character_alphabet if item != '-15']
    site_patterns = list(product(character_alphabet,repeat = N_taxa))
    site_patterns = list(map('_'.join,site_patterns))
    print(f"\nTotal site patterns: {len(site_patterns)}")
    dtype = [tuple([i,"f8"]) for i in site_patterns]
    struc_array = np.zeros(N_alns, dtype=dtype)
    #dic_patterns = dict.fromkeys(site_patterns,np.repeat(0,N_alns))
    for i in range(N_alns):
        aln = array_in[i,:,:]
        aln = aln[:,(aln != -15).any(axis=0)]
        Aln_length = aln.shape[1] 
        for s in range(Aln_length):     
            array_key = '_'.join((map(str,aln[:,s].flatten())))
            struc_array[array_key][i]+=1/Aln_length  
    return(struc_array.view((float,len(struc_array.dtype.names))))      


#Regression: inference of branch lengths of unrooted trees from site patterns  
def build_MLP_brl(X_train,Y_train,droput_rates,batch_sizes,X_test):
    
    # Length of feature vector, i.e. number of patterns 
    N_patterns=X_train.shape[1]
    
    #Number of branches of a tree
    N_branch = Y_train.shape[1]
    
    visible_layer = Input(shape=(N_patterns,))
    hidden1 = Dense(1000,activation='relu')(visible_layer)
    drop1 = Dropout(rate=droput_rates)(hidden1)
    hidden2 = Dense(1000,activation='relu')(drop1)
    drop2 = Dropout(rate=droput_rates)(hidden2)
    hidden3 = Dense(1000,activation='relu')(drop2)
    drop3 = Dropout(rate=droput_rates)(hidden3)
    output = Dense(N_branch, activation='linear')(drop3)
    
    #Randomly pick 100 from test dataset 
    trace_data = X_test[np.random.choice(range(200),100,replace=False)]
    #Trace model helper
    model_trace = Model(inputs=visible_layer, outputs=[drop3])
    trace = m_phate.train.TraceHistory(trace_data, model_trace)
    #Save loss values
    history = keras.callbacks.History()
    #Model compilation
    model_mlp = Model(inputs=visible_layer, outputs=output)
    model_mlp.compile(loss='mean_squared_error',optimizer='adam',metrics=['mae','mse'])
    
    
    #Print model
    print(model_mlp.summary()) 
    
    model_mlp.fit(x=X_train,y=Y_train,batch_size=batch_sizes,callbacks=[trace,history],epochs=200,verbose=1,shuffle=True,validation_split=0.1)
    return model_mlp, trace, history


def main():
    parser = argparse.ArgumentParser(description='Keras run')
    parser.add_argument( '--tr', help = "Training MSAs dataset in npy",nargs='+',dest='TRAIN')
    parser.add_argument( '--te', help = "Test MSAs dataset in npy",nargs='+',dest='TEST')
    parser.add_argument( '--trl', help = "Training branch lengths in csv",nargs='+',dest='LABTRAIN')
    parser.add_argument( '--tel', help = "Test branch lengths in csv",nargs='+',dest='LABTEST')
    parser.add_argument( '--dr', help = "Dropout rate e.g. 0.15",dest='DROP',type=float,default=0.15)
    
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
    model_mlp_reg, my_trace, my_history = build_MLP_brl(X_train=X_train,Y_train=Y_train,droput_rates=args.DROP,batch_sizes=100, X_test=X_test)
    #np.savetxt("randomseed.mlp.mphate.dropout_{}.txt".format(args.DROP),np.random.get_state()[1][0],fmt='%f')
    print("Running M-PHATE")
    #Extract trace data
    trace_data = np.array(my_trace.trace)
    #Apply M-PHATE
    #Extract epoch ID and N neurons 
    epoch = np.repeat(np.arange(trace_data.shape[0]), trace_data.shape[1])
    n_neurons = trace_data.shape[1]
    #The train loss for each element of the flattened trace
    loss = np.repeat(my_history.history['mse'], n_neurons)
    np.savetxt("loss_train.mlp.mphate.dropout_{}.txt".format(args.DROP),loss,fmt='%f')
    #The validation loss for each element of the flattened trace
    val_loss = np.repeat(my_history.history['val_mse'], n_neurons)
    np.savetxt("loss_val.mlp.mphate.dropout_{}.txt".format(args.DROP),val_loss,fmt='%f')
    #Visualization parameters 
    param_name = 'intraslice_knn'
    params = [1, 2, 5, 10]
    
    for param in params:
        m_phate_op = m_phate.M_PHATE(verbose=1, **{param_name : param})
        m_phate_data = m_phate_op.fit_transform(trace_data)
        fig = scprep.plot.scatter2d(m_phate_data,c=epoch,ticks=True,title="Epoch",label_prefix="M-PHATE")
        fig.figure.savefig("mlp_{}_{}_dropout_{}_epoch.png".format(param_name, param, args.DROP),dpi=300)
        fig = scprep.plot.scatter2d(m_phate_data, c=loss, ticks=True, title='Training loss', label_prefix="M-PHATE",cmap_scale="sqrt")
        fig.figure.savefig("mlp_{}_{}_dropout_{}_loss.png".format(param_name, param, args.DROP),dpi=300)
        fig = scprep.plot.scatter2d(m_phate_data, c=val_loss, ticks=True, title='Validation loss', label_prefix="M-PHATE",cmap_scale="sqrt")
        fig.figure.savefig("mlp_{}_{}_dropout_{}_valloss.png".format(param_name, param, args.DROP),dpi=300)
        m_phate_op.set_params(n_components=3)
        m_phate_data = m_phate_op.transform()
        np.savetxt("dims.mlp.mphate.dropout_{}.txt".format(args.DROP),m_phate_data,fmt='%f')
        scprep.plot.rotate_scatter3d(m_phate_data, c=epoch,
                             title='Epoch', figsize=(10,8),
                             ticks=False, label_prefix="M-PHATE",
                             legend_anchor=(1.05, 1),
                             filename="mlp_{}_{}_dropout_{}_epoch.gif".format(param_name, param, args.DROP), fps=25, rotation_speed=60,
                             )
        scprep.plot.rotate_scatter3d(m_phate_data, c=loss,
                             title='Training loss', figsize=(10,8),
                             ticks=False, label_prefix="M-PHATE",
                             legend_anchor=(1.05, 1),
                             filename="mlp_{}_{}_dropout_{}_loss.gif".format(param_name, param, args.DROP), fps=25, rotation_speed=60,
                             )
        scprep.plot.rotate_scatter3d(m_phate_data, c=val_loss,
                             title='Validation loss', figsize=(10,8),
                             ticks=False, label_prefix="M-PHATE",
                             legend_anchor=(1.05, 1),
                             filename="mlp_{}_{}_dropout_{}_valloss.gif".format(param_name, param, args.DROP), fps=25, rotation_speed=60,
                             )
    
    #Save M-PHATE results 
    #np.save("phate.mlp.trace",trace_data)
    #np.save("phate.mlp.result",m_phate_data)
    #Evaluate model
    #print("Evaluating with best weights")
    #evals_reg = model_mlp_reg.evaluate(X_test,Y_test,batch_size=100, verbose=1, steps=None)
    #bls = model_mlp_reg.predict(X_test,batch_size=100, verbose=1, steps=None)
    #np.savetxt("brls.evaluated.mlp.txt",evals_reg,fmt='%f')
    #np.savetxt("brls.predicted.mlp.txt",np.exp(bls),fmt='%f')
    #Saving model
    #print("\nSaving keras trained model")
    #model_mlp_reg.save("keras_model_mphate_mlp.h5")
    tf.keras.utils.plot_model(model_mlp_reg, to_file='model_mphate_mlp.png', show_shapes=True)
       
    
if __name__ == "__main__":
    main() 