'''
Created by Ryan Schenck
20 March 2018

This script is designed to work on the new dataset of final variant calls from the PCAWG dataset for predictions.
'''

try:
    from optparse import OptionParser
    import logging
    import pickle
    from sklearn.metrics import roc_auc_score, roc_curve, auc
    import datetime
    from scipy import interp
    from functools import wraps
    import time
    import sys
    import os
    import glob
    from collections import OrderedDict

    import keras as ks
    import numpy as np
    import h5py
except Exception as e:
    print(e, file=sys.stdout)
    sys.exit("ERROR: Unable to load dependencies.")

'''
Model To Be Used for Evaluations!!!
Local: /Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/
Cluster: /well/wedge/rschenck/DPhilRotation1/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/
'''

def OptionParsing():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-m', '--model', dest="InputDir", default="/Users/schencro/Desktop/Oxford/Rotation_1/CNN/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59/", help="Directory where the model is located to make the predictions.")
    parser.add_option('-o', '--output', dest="output", default="/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/WTPreds.txt", help="Output directory")
    parser.add_option('-i', '--input', dest='infasta', default='/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/WTseqs.fasta', help="Input fasta file")
    parser.add_option('-d', '--debug', dest='debug', default=False, action='store_true', help="Run to do small batch predictions.")
    (options, args) = parser.parse_args()
    return (options, parser)

def fn_timer(function):
    '''
    Use this as a wrapper at the top of any function you want to get run time information about.

    :param function: Function of interest.
    :return: A function to wrap around a function.
    '''
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("INFO: Total time running %s: %s minutes" %
               (function.__name__, str(round((t1-t0)/60.,2)))
               )
        return result
    return function_timer

def LoadModel(Options):
    model_to_load = glob.glob(Options.InputDir + "*.modelConfig.yaml")[0] # Should only be one
    weights_to_load = glob.glob(Options.InputDir + "*.modelWeights.h5")[0]

    # Load Model Configuration from yaml
    with open(model_to_load, 'r') as inModel:
        loadedModel = inModel.read()

    loaded_model = ks.models.model_from_yaml(loadedModel)
    print("Model Configuration Loaded.", file=sys.stdout)

    # Load Weights from h5
    loaded_model.load_weights(weights_to_load)
    print("Model weights loaded.", file=sys.stdout)

    return(loaded_model)

def dna_one_hot(seq, seq_len=None, flatten=True):
    '''
    :param seq: Input sequence read from fasta file.
    :param seq_len: Length of longest sequence
    :param flatten: Whether or not to flatten into a column vector (May Change)
    :return: seq_vec: A Flattened column vector (May change)
    '''
    if seq_len == None:
        seq_len = len(seq)
        seq_start = 0
    else:
        if seq_len <= len(seq):
            # trim the sequence
            seq_trim = (len(seq)-seq_len) // 2
            seq = seq[seq_trim:seq_trim+seq_len]
            seq_start = 0
        else:
            seq_start = (seq_len-len(seq)) // 2

    seq = seq.upper()

    seq = seq.replace('A','0')
    seq = seq.replace('C','1')
    seq = seq.replace('G','2')
    seq = seq.replace('T','3')

    # map nt's to a matrix 4 x len(seq) of 0's and 1's.
    #  dtype='int8' fails for N's
    seq_code = np.zeros((seq_len,4), dtype='float16')
    for i in range(seq_len):
        if i < seq_start:
            seq_code[i:] = 0.25
        else:
            try:
                seq_code[i,int(seq[i-seq_start])] = 1
            except:
                seq_code[i:] = 0.25

    # flatten and make a column vector 1 x len(seq)
    seq_vec = seq_code.flatten()[None,:]
    # seq_vec = seq_code
    # print(seq_vec.reshape(600,4)) # This is the appropriate shape after flattening...

    return seq_vec

def RunPredictions(seq_vecs, model, Options):
    # Step 1 separate headers and chunks
    seqs = []
    heads = []
    for head in seq_vecs:
        seqs.append(seq_vecs[head])
        heads.append(head)

    seqStack = np.vstack(seqs)
    seqsReady = seqStack.reshape(seqStack.shape[0], 600, 4)

    preds = model.predict(seqsReady)
    #
    # # outLine = (heads[0] + '\t' + '\t'.join([repr(i) for i in preds.tolist()[0]]) +'\n')
    outLine = heads[0] + '\t' + repr(np.mean(preds)) + '\t' + repr(np.std(preds)) + '\n'
    with open(Options.output, 'a') as outFile:
        outFile.write(outLine)


@fn_timer
def ProcessFasta(Options, model):
    #### Takes two and a half hours to process one full PCAWG fasta set (2.25 * 2 for both WT and Mut) with chunk size of 250000
    i = 0
    seq = ''
    chunkSize = 1
    seq_vecs = OrderedDict()
    with open(Options.infasta, 'r') as inFile:
        for line in inFile:
            if line[0] == '>':
                if seq:
                    i += 1
                    # Main code chunck
                    ###
                    if len(seq_vecs)==chunkSize-1:
                        seq_vecs[header] = dna_one_hot(seq, 600) # Not flattened
                        # do something
                        RunPredictions(seq_vecs, model, Options)
                        # if i==1:
                        #     break
                        # TODO Make predictions once there are 250000 seqs put together
                        ##
                        seq_vecs=OrderedDict()
                    else:
                        seq_vecs[header] = dna_one_hot(seq, 600)

                    ###
                header = line[1:].rstrip()
                seq = ''
            else:
                seq += line.rstrip()

@fn_timer
def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('DataExtractionForMutagenesis.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()
    os.system('touch %s' % (Options.output))
    model = LoadModel(Options)

    # model.summary()
    ProcessFasta(Options, model)

if __name__=='__main__':
    main()