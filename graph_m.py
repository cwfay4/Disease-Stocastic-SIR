import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import itertools as tools

iterations = []
susc = []
infe = []
imm = []
dead = []
susc_err = []
infe_err = []
imm_err = []
dead_err = []

st=input("Enter start: ")
start=int(st)
stop=start+500

with open('Dconv.csv', 'r') as csv_file:  #importing the csv file
    csv_reader = csv.reader(csv_file)  #reading the file
    
    for line in tools.islice(csv_file,start,stop):
        currentline = line.split(",") #splits the lines into pieces based on commas so we can sort it
        x=currentline[0]  #this could work {iterations}
        iterations.append(float(x))
        y1=currentline[1] #{susceptible}
        susc.append(float(y1))
        y1_err=currentline[2] #{error on our susceptible average}
        susc_err.append(float(y1_err))
        y3=currentline[3] #{infected}
        infe.append(float(y3))
        y3_err=currentline[4]
        infe_err.append(float(y3_err)) #{error on our infected average}
        y5=currentline[5] #{immune}
        imm.append(float(y5))
        y5_err=currentline[6]
        imm_err.append(float(y5_err)) #{error on our immune average}
        y7=currentline[7] #{dead}
        dead.append(float(y7))
        y7_err=currentline[8]
        dead_err.append(float(y7_err))
        print(x,y1,y1_err,y3,y3_err,y5,y5_err,y7,y7_err)

imm_i = [float(i) for i in imm]

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']

plt.errorbar(iterations,susc,yerr=susc_err,label='Susceptible',ecolor='r',elinewidth=0.5,barsabove=True)
plt.errorbar(iterations,infe,yerr=infe_err,color='darkorange',label='Infected',ecolor='mediumpurple',elinewidth=0.5,barsabove=True)
if sum(imm_i) != 0:
    plt.errorbar(iterations,imm,yerr=imm_err,color='g',label='Immune',ecolor='maroon',elinewidth=0.5,barsabove=True)
plt.errorbar(iterations,dead,yerr=dead_err,color='black',label='Dead',ecolor='yellow',elinewidth=0.5,barsabove=True)

plt.title('100000 2 .001 .7 .2 0 .1',fontsize=14)
plt.xlabel('Iterations',fontsize=12)
plt.ylabel('Fraction of Connected Population',fontsize=12)

path='C:\PHYS3500\CompPhys\\'
plt.legend(loc='upper right',fontsize='large')
plt.savefig(os.path.join(path,'-sclf_100k_1.2_.001_.05_.1_.05.pdf'),dpi=500)
plt.show()