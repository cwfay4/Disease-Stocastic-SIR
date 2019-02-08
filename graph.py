import numpy as np
import matplotlib.pyplot as plt
import os
import csv

iterations = []
susc = []
infe = []
imm = []
dead = []

with open('Disease_Conv.csv', 'r') as csv_file:  #importing the csv file
    csv_reader = csv.reader(csv_file)  #reading the file
    
    for line in csv_file: #reading the csv file, sorting it so we can graph it
        currentline = line.split(",") #splits the lines into pieces based on commas so we can sort it
        x=currentline[0]  #this could work {iterations}
        iterations.append(x)
        y1=currentline[1] #{susceptible}
        susc.append(y1)
        y2=currentline[2] #{infected}
        infe.append(y2)
        y3=currentline[3] #{immune}
        imm.append(y3)
        y4=currentline[4] #{dead}
        dead.append(y4)
        #print(x,y1,y2,y3,y4)

imm_i = [float(i) for i in imm]

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']

plt.plot(iterations,susc,label='Susceptible')
plt.plot(iterations,infe,label='Infected')
if sum(imm_i) != 0:
    plt.plot(iterations,imm,label='Immune')
plt.plot(iterations,dead,label='Dead',color='r')

plt.title('-hcm 11 4.2 100000 1.5 .001 .8 .3 0 .1',fontsize=14)
plt.xlabel('Iterations',fontsize=12)
plt.ylabel('Fraction of Connected Population',fontsize=12)

plt.legend(loc='upper right',fontsize='large')
#plt.savefig(os.path.join(path,'hcm_trial.pdf'),dpi=300)
plt.show()