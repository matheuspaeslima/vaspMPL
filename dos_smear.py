#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import os
import sys

print("###############################################")
print("#                                             #")
print("#  Title : dos_smear.py                       #")
print("#  Author: Prof. Matheus Paes Lima            #")
print("#  Date  : 9Jul23                             #")
print("#  Descr.: This code smears density-of-states #")
print("#          (DOS) data, originally calculated  #")
print("#          using VASP code and processed by   #") 
print("#          vaspkit.                           #")
print("#                                             #")
print("###############################################")

#============= input ===============================#
fname=input('Provide the input file name:')

#Output file name
fname_split=fname.rsplit('.',1)
fname_out=fname_split[0]+'_smear'+'.'+fname_split[1]
if os.path.exists(fname_out):
	print(   "ERROR: The file ",fname_out," exists")
	sys.exit("===== stopping the code =====")

smear=float(input('Provide the smearing parameter in eV:'))
if (smear<0.01):
	print(   "ERROR: This program requires smearing larger than 0.01eV")
	sys.exit("===== stopping the code =====")
if (smear>0.5):
	print(   "WARNING: This program suggests smearing larger than 0.5eV")
	
#============= Reading Data Frame ==================#
df = pd.read_csv(fname, delimiter=f'\s+')

#== smearing the dataframe (data will be replaced) =#
cols=df.columns
x = df[cols[0]]
h=x[1]-x[0]
sigma=smear/h
for col in cols:
    y = df[col]
    y_smear = gaussian_filter(y, sigma)
    df[col]=y_smear

#============= Writing data to file ================#

#test: df.to_csv(fname_out,sep=' ',index=False,float_format='%.5f')

# Write the formatted string to a file
with open(fname_out, 'w') as file:
	#header === 
	for item in cols:
		formatted_item = '{:>13}'.format(item)
		file.write(formatted_item)
	file.write('\n')

	#data ===
	# Convert dataframe to a string with desired formatting
	data_string = df.to_string(index=False, header=False, float_format=lambda x: f"{x:12.5f}")
	file.write(data_string)
