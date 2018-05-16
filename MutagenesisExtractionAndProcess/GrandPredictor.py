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
    parser.add_option('-o', '--output', dest="output", default="Predictions/", help="Output directory")
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

@fn_timer
def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('DataExtractionForMutagenesis.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()

    model = LoadModel(Options)

    # model.summary()


if __name__=='__main__':
    main()