####################################################################################
# this program use tensorflow to learn the reconstruction of sparse odor mixture
# comparison of reconstruction with different basal activity when considering spontaneous
# activity and odor-evoked inhibition
# The first layer weightmatrix is loaded from previous simulations
# USAGE:
# RUN in terminal: python olfRconInhi N M noise L r0
# N       number of odorant, default 50
# M       number of receptors, default 10
# noise   standard deviation of response noise, 0.05
# L       number of hidden layers, an integer
# r0      relative basal activity, 0 ~ 1

# version 0.1
# Initially writen by Qianyi Li, modified by Shanshan Qin, 12/29/2018
######################################################################################
from os.path import dirname, join as pjoin
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
    N = int(sys.argv[1])
else:
    N = 50   #default

if total > 2:
    M = int(sys.argv[2])
else:
    M = 10 #default

if total > 3:
    noise = float(sys.argv[3])
else:
    noise = 0.01 #default

if total > 4:
    L = int(sys.argv[4])
else:
    L = 2

if total > 5:
    r0 = float(sys.argv[5])
else:
    r0 = 0.95
    

# Training Parameters
learning_rate = 0.01
num_steps = 10000
batch_size = 1000
thr = 1  # this is used in error function, only when considering ln(thr + c)
display_step = 100

# set up network parameters
num_input = N   # number of odoratns
num_output = M    # number of receptors
num_hidden = 100  
sample_size = 1000000
sparsity = 2
sigma = 2

# load the matrix data
dFolder = '../data/inhi_N50M10S2sig2_1013'
fName = 'gcmiInhi_N50_R10_S2_sig2_alp' + str(r0) + '_frac_2018-10-13.mat'
dPath = pjoin(dFolder,fName)
mat = scio.loadmat(dPath)
allW = mat["allMat"]
allSign = mat["allSign"]         #all the signs of matrix elements
allfmin = mat["allfmin"]

# define network shape
netShape = np.ones(L+3,dtype=int)*num_output
netShape[1:(L+1)] = np.ones(L, dtype = int)*num_hidden
netShape[-2:] = [num_input,num_input]
# 
repeat = len(allfmin);  # repeats of one basal activity
testLost = np.zeros([repeat])
trainLost = np.zeros([repeat])

X = tf.placeholder("float", [None, num_input])
X_pr = tf.placeholder("float", [None, num_input])
odor_binary_tf = tf.placeholder("float", [None, num_input])

#%%
#def weight_variable(shape):
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


def decoder(x):
    activation = x
    #zs = []
    for w, b in zip(weights[:-2], bias[:-2]):
        activation = tf.nn.sigmoid(tf.add(tf.matmul(activation,w), b))
    reconstruct = tf.add(tf.matmul(activation, weights[-2]), bias[-2])
#    classification = tf.add(tf.matmul(reconstruct, weights[-1]), bias[-1])
    classification = tf.add(tf.multiply(reconstruct, weights[-1]), bias[-1])
    return reconstruct, classification

def resp(x,we,wi,r0):
    alpha = (1-r0)/r0
    layer_1 = 1/(1+alpha*(1+tf.matmul(x,wi))/(1 + tf.matmul(x,we)))
    layer_1 = tf.add(layer_1, noise*tf.random_normal([num_output]))
    return layer_1
     
#%%
input_data, input_nonzero_data, input_binary = genTrainData(num_input, num_output, sample_size, sigma)

for k in range(0, repeat):
    # store the loss function list
        loss_list = []
        # reshape the excitatory and inhibitory matrix
        W = allW[:,k].reshape(num_input,num_output)
        S = allSign[:,k].reshape(num_input,num_output)
        we = np.zeros([num_input,num_output])
        wi = np.zeros([num_input,num_output])
        we[S==1] = W[S==1]
        wi[S==-1] = W[S==-1]
        we = tf.cast(we,tf.float32)
        wi = tf.cast(wi,tf.float32)
        
        encoder_op = resp(X,we,wi,r0)                            
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
    
    # save some data for future comparison
save_loss = './loss_inhi_' +  str(N) + 'M' + str(M) + 'noise' + str(noise) + '_sp' + str(r0) + '_L' + str(L) + '.mat'
scio.savemat(save_loss, {'testLost':testLost,'trainLost':trainLost, 'prediction':prediction, 'original':original})
    # plot the data
#    plt.subplot(1, 2, 1)
#    plt.imshow(prediction, vmin=-5, vmax=5)
#    plt.subplot(1, 2, 2)
#    plt.imshow(original, vmin=-5, vmax=5)
#    plt.show()
#    plt.subplot(1, 2, 1)
#    plt.hist(weight_pr, bins=20)
#    plt.subplot(1, 2, 2)
#    plt.imshow(sess.run(weights['encoder_h1']), vmin=-3, vmax=3)
#    plt.show()
#    
#    plt.plot(loss_list)