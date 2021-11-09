import os
import numpy as np

polymers =[]
timestep = 1000
atoms = 400
maxtime = 4001000
l = 100
inputs= ["dump.DNAunzipping1.3.kamograph", "dump.DNAunzipping1.4.kamograph",
"dump.DNAunzipping1.5.kamograph", "dump.DNAunzipping1.8.kamograph"]


averagekamos = open(input + ".average.txt", "w")
for input in inputs:
    total = 0
    with open("/home/s1718990/DNAUnzipping/outputs/"+input+".txt") as fp:
        while timestep <= maxtime:
            print (timestep)
            fp.seek(0)
            for i, line in enumerate(fp):
                if (i >= (timestep/1000 * atoms) and i < ((timestep/1000 + 1) * atoms)):

                    newline = line.split()
                    total += newline[2]

        total = float(total)/float(atoms)
        averagekamos.write(str(timestep) + " " + str(total)+ '\n')

            print(timestep)
            timestep += 1000
