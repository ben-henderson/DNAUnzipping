import math
import os
import numpy as np

polymers =[]
results = []
timestep = 0
atoms = 400
maxtime = 4001000
lx = 100
input = "dump.DNAunzipping1.0"
polymers = np.zeros((atoms, 3),dtype=float)
with open("/home/ben/github/DNAUnzipping/outputs/"+input) as fp:
    while timestep <= maxtime:
        for i, line in enumerate(fp):
            if i >= (timestep*(atoms+9))+9 and i < (timestep+1)*(atoms+9):
                newline = line.split()
                #lines.append(newline)
                atomid = int(newline[0])-1
                x = float(newline[2])
                y = float(newline[3])
                z = float(newline[4])
                polymers[atomid, 0] = x*100
                polymers[atomid, 1] = y*100
                polymers[atomid, 2] = z*100


        for x in range(int((atoms/2)-2)):
            if (x == 5): print(timestep, polymers[x,1], polymers[x+200, 1])
            dx = math.abs(polymers[x,0] - polymers[x+200,0])
            dy = math.abs(polymers[x,1] - polymers[x+200,1])
            dz = math.abs(polymers[x,2] - polymers[x+200,2])

            #angle = math.atan2(dy/dx)
            #if angle < 2.5:
                #result = [timestep*1000, lines[x*2][1], 0]
#
#            else:
#                result = [timestep*1000, lines[x*2][1], 1]
#
#            results.append(result)
#            timestep =+ 1000


#print (results)
