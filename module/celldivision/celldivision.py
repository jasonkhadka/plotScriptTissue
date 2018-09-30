####################################################################################################################################
#                     MODULE TO DO CELLDIVISION
####################################################################################################################################
#importing all the libraries
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
#pd.set_option('precision',10)
#pd.set_option("display.max_columns",200)
#pd.set_option("display.max_rows",2000)
import sys 
sys.path.append("/home/jkhadka/plantdev")
sys.path.append("/home/jkhadka/plantdev/python_quadedge")
import quadedge as qd
import Quadedge_lattice_development as latdev
import centered_lattice_generator as latgen
import matplotlib.pyplot as plt
#from IPython.display import display, HTML
import nlopt #non linear optimizer
import argparse #argument parser, handles the arguments passed by command line
import os # to make directory
import time
#importing functions that the simulations rely on -- make sure to copy simulation_functions.py to the same directory as this file
from simulation_functions import *
from datetime import datetime
##################################################################################################################################################
programstart = time.time()
print "This is a test for Growth !!"
algorithm = "SBPLX"
print "The Used Algorithm :: ", algorithm