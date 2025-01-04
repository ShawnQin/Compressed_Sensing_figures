##########################################################################
# This program use tensorflow to learn the reconstruction of sparse odor mixture
# We want to demonstrate that an optimal matrix that achieves maximum entropy coding
# performance. We change the coding scheme by varying the sparisty of sensitivity matrix
# compared with learning of both encoding and decoding

# RUN in terminal: python olfRconInhi N M noise L r0
# N       number of odorant, default 50
# M       number of receptors, default 10
# noise   standard deviation of response noise, 0.05
# L       number of hidden layers, an integer
# firstW  plastic or hard wired
#
# first written by Qianyi Li
# version 0.2
# the last layer, the classification layer has been revised
# modified by Shanshan Qin
# last revision on 12/29/2018
######################################################################################

import tensorflow as tf
import numpy as np
import random
#import matplotlib.pyplot as plt
import scipy.io as scio
import sys

# parse the input argument
total = len(sys.argv)
cmdargs = str(sys.argv)

if total > 1:
    ix = int(sys.argv[1])
else:
     ix = 17
     
if total > 2:
    N = int(sys.argv[2])
else:
    N = 20   #default

if total > 3:
    M = int(sys.argv[3])
else:
    M = 9 #default

if total > 4:
    noise = float(sys.argv[4])
else:
    noise = 0.01 #default

if total > 5:
    L = int(sys.argv[5])
else:
    L = 2

if total > 6:
    firstW = sys.argv[6]
else:
    firstW = "hardWire"
    

# Training Parameters
learning_rate = 0.01
num_steps = 10000
batch_size = 1000
thr = 1  # this is used in error function where ln(thr + c) is compareds
display_step = 100

# set up network parameters
num_input = N   # number of odoratns
num_output = M    # number of receptors
num_hidden = 100  
sample_size = 1000000
sparsity = 2
sigma = 2

# define network shape
netShape = np.ones(L+3,dtype=int)*num_output
netShape[1:(L+1)] = np.ones(L, dtype = int)*num_hidden
netShape[-2:] = [num_input,num_input]
#sp2 = 0.45      # sparsity of first layer weight


repeat = 1;
allSp = np.arange(0.05,1.05,0.05)   # sparsity of sensitivity matrix
testLost = np.zeros([repeat])
trainLost = np.zeros([repeat])
sp2 = allSp[ix]  #select the sp

X = tf.placeholder("float", [None, num_input])
X_pr = tf.placeholder("float", [None, num_input])
odor_binary_tf = tf.placeholder("float", [None, num_input])

# weight matrix at each layer of the decoder
def weight_variable(shape):
    w =[]
    for m,n in zip(shape[:-2],shape[1:-1]):
        w.append(tf.Variable(tf.random_normal([m,n],stddev=0.1)))
    w.append(tf.Variable(tf.random_normal([1,shape[-1]],stddev=0.1)))
    return w
def bias_variable(shape):
    return [tf.Variable(tf.random_normal([1,m], stddev=0.1)) for m in shape[1:]]

weights = weight_variable(netShape)
bias = bias_variable(netShape)

# if first layer is also platic
if firstW == "plastic":
    W1 = tf.Variable(tf.random_normal([N,M],stddev=0.1))

# random input matrix with sparsity
def sparseWeight(num_input, num_output, sp2, sigma):
    weight = np.exp(np.random.multivariate_normal(-np.ones(num_output), sigma**2*np.eye(num_output), num_input))
    zero_ind = random.sample(list(np.arange(0, num_input*num_output)), int(round(sp2*num_input*num_output)))
    weight = np.reshape(weight, [num_output*num_input, 1])
    weight[zero_ind] = 0
    weight = np.reshape(weight, [num_input, num_output])
    weight = np.float32(weight)
    return weight

# generating training data
def genTrainData(num_input, num_output, sample_size, sigma):
    input_data = np.zeros([num_input, sample_size])
    input_nonzero_data = np.zeros([num_input, sample_size])
    input_nonzero = np.exp(np.random.multivariate_normal(np.zeros(sparsity), sigma**2*np.eye(sparsity), sample_size))
    for j0 in range(0, sample_size):
        a = random.sample(list(np.arange(0, num_input)), sparsity)
        input_data[a, j0] = input_nonzero[j0, :]
        input_nonzero_data[a, j0] = np.log(input_nonzero[j0, :])
        
    input_data = input_data.transpose()
    input_nonzero_data = input_nonzero_data.transpose()
    input_binary = input_data > 0
    return input_data, input_nonzero_data, input_binary

# the encoder, first hidden layer, with response noise
def encoder(x):
    # Encoder Hidden layer with sigmoid activation #1
#    layer_1 = tf.matmul(x, tf.exp(weights['encoder_h1']))/(1+tf.matmul(x, tf.exp(weights['encoder_h1'])))
    layer_1 = tf.matmul(x, tf.exp(W1))/(1+tf.matmul(x, tf.exp(W1)))
    layer_1 = tf.add(layer_1, noise*tf.random_normal([num_output]))
    return layer_1

def decoder(x):
    activation = x
    #zs = []
    for w, b in zip(weights[:-2], bias[:-2]):
        activation = tf.nn.sigmoid(tf.add(tf.matmul(activation,w), b))
        #zs.append(z)       
    reconstruct = tf.add(tf.matmul(activation, weights[-2]), bias[-2])
#    classification = tf.add(tf.matmul(reconstruct, weights[-1]), bias[-1])
    classification = tf.add(tf.multiply(reconstruct, weights[-1]), bias[-1])
    return reconstruct, classification

# used for fixed input weight, non-plastic
def resp(x,weight):
    layer_1 = tf.matmul(x, weight)/(1+tf.matmul(x, weight))
    layer_1 = tf.add(layer_1, noise*tf.random_normal([num_output]))
    return layer_1

#%%
input_data, input_nonzero_data, input_binary = genTrainData(num_input, num_output, sample_size, sigma)

#for s in range(0, len(allSp)):
for k in range(0, repeat):
    # store the loss function list
        loss_list = []
        #sp2 = allSp[s]
        if firstW == "hardWire":
            print("sparisty of w is: %f" % sp2)
            weight = sparseWeight(num_input, num_output, sp2, sigma)
            encoder_op = resp(X,weight) 
        elif firstW == "plastic":
            encoder_op = encoder(X)  # for plastic input plasticity
        else:
            print("the first layer type is wrong!")    
        
#        encoder_op = resp(X,weight)                # for fixed input plasticity
                            
        decoder_op = decoder(encoder_op)[0]
        class_op = decoder(encoder_op)[1]

        y_pred = decoder_op
        y_true = X_pr
        y_class = class_op
        pred = tf.multiply(y_pred, odor_binary_tf)

        # fixed first layer and changing sparsity
        #class_err = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=odor_binary_tf, logits=class_op))
        class_err = tf.reduce_mean(tf.maximum(class_op,0) - tf.multiply(class_op, odor_binary_tf) \
                    + tf.log(1 + tf.exp(-tf.abs(class_op))))
        recons_err = tf.reduce_mean(tf.pow(tf.multiply(y_pred, odor_binary_tf) -
                            tf.multiply(y_true, odor_binary_tf), 2))\
                            * num_input / sparsity
        loss = recons_err + class_err
        optimizer = tf.train.AdamOptimizer(learning_rate).minimize(loss)
        
        #initialize
        init = tf.initialize_all_variables()
        
        with tf.Session() as sess:
    
        # Run the initializer
            sess.run(init)
            
            # Training
            for i in range(1, num_steps+1):            
                # Get the next batch of input data
                j = divmod(i, int(round(sample_size/batch_size)))[1]+1
                batch_x = input_data[batch_size*(j-1):batch_size*j, 0:]
                binary_x = input_binary[batch_size*(j-1):batch_size*j, 0:]
                non_zero = input_nonzero_data[batch_size*(j-1):batch_size*j, 0:]
                 #print(batch_x)
        
                # Run optimization op (backprop) and cost op (to get loss value)
                _, l, p, t = sess.run([optimizer, loss, recons_err, class_err], feed_dict={X: batch_x, odor_binary_tf: binary_x, X_pr: non_zero})
                loss_list.append(l)  # store the loss function
                # Display logs per step
                if i % display_step == 0 or i == 1:
                    print('Step %i: Minibatch Loss: %f Reconstruction Loss: %f Classification Loss: %f' % (i, l, p, t))
                    #weight_pr = np.reshape(sess.run(weights['encoder_h1']), [1, num_output*num_input]).transpose()
            
            trainLost[k] = np.mean(loss_list[-100:])
            
            #generating testing data
            test_data, test_nonzero_data, test_binary = genTrainData(num_input, num_output, 10000, sigma)       
            testLost[k] = sess.run(loss,feed_dict={X: test_data, odor_binary_tf: test_binary, X_pr: test_nonzero_data})
            
            # small set used to make scatter plots
            test_set = test_data[0:1000, 0:]
            test_set_true = test_nonzero_data[0:1000, 0:]
            test_set_bin = test_binary[0:1000, 0:]
            prediction = sess.run(pred, feed_dict={X: test_set, odor_binary_tf: test_set_bin, X_pr: test_set_true})
            original = sess.run(y_true, feed_dict={X: test_set, X_pr: test_set_true, odor_binary_tf: test_set_bin})
    
    # save some data into mat file and will be used for future comparison
save_loss = './loss_' + firstW + '_'+  str(N) + 'M' + str(M) + 'noise' + str(noise) + '_sp' + str(sp2) + '_L' + str(L) + '.mat'
scio.savemat(save_loss, {'testLost':testLost,'trainLost':trainLost, 'prediction':prediction, 'original':original})