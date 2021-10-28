import os
import numpy as np

polymers =[]
timestep = 0
atoms = 400
maxtime = 4001000
l = 100
input = "dump.DNAunzipping2.5"
polymers = np.zeros((atoms, 3),dtype=float)
kamo = open(input + ".kamograph.txt", "w")
with open("/home/ben/github/DNAUnzipping/outputs/"+input) as fp:
    while timestep <= maxtime:
        for i, line in enumerate(fp):
            if (i >= ((timestep/1000 * atoms) + 9) and i < ((timestep/1000 + 1) * atoms)):
            #if line in lines_to_read:
                newline = line.split()
                #lines.append(newline)
                atomid = int(newline[0])-1
                x = float(newline[2])
                y = float(newline[3])
                z = float(newline[4])
                polymers[atomid, 0] = x*100
                polymers[atomid, 1] = y*100
                polymers[atomid, 2] = z*100
                #print(x, y, z)



        for x in range(int((atoms/2))):
            #if (x == 5): print(timestep, x, polymers[x+200, 1])
            dx = np.abs(polymers[x + 200,0] - polymers[x,0])
            dy = np.abs(polymers[x + 200,1] - polymers[x,1])
            dz = np.abs(polymers[x + 200,2] - polymers[x,2])
            #print (dx)
            if (dx > l/2):
                dx = l - dx

            if (dy > l/2):
                dy = l - dx

            if (dz > l/2):
                dz = l - dz

            dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if (dist < 2.5):
                result = 0
            else:
                result = 1
            #print (dx, dx ** 2)
            #print (dy, dy ** 2)
            #print (dist)
            kamo.write(str(timestep) + " " + str(x) + " " + str(result)+ '\n')
        #print(timestep)
        timestep += 1000
