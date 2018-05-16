'''
Created by Ryan Schenck
20 March 2018

This script is designed to work on the new dataset of final variant calls from the PCAWG dataset.
'''
import sys
import os
import glob
from optparse import OptionParser
from collections import OrderedDict
import time
import subprocess
from functools import wraps
import numpy as np
import gzip
import pickle
import subprocess

'''
Model To Be Used for Evaluations!!!
/well/wedge/rschenck/DPhilRotation1/Model/AllENCODEnocancer.23Feb2018.batch64.2018-02-23.21.59
'''

def OptionParsing():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-m', '--model', dest="model", default=None, help="Directory where the model is located to make the predictions.")
    parser.add_option('-o', '--output', dest="output", help="Directory where the vcf files are stored. This will process all files in a directory.")
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


@fn_timer
def main():
    FilePath = os.path.dirname(os.path.abspath(__file__))
    localpath = os.path.abspath(__file__).rstrip('DataExtractionForMutagenesis.py')  # path to scripts working directory

    (Options, Parser) = OptionParsing()

    pass

if __name__=='__main__':
    main()