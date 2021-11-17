import os
import numpy as np

atoms = 400
maxtime = 4001000
eqtime = 2000000
inputs= ["dump.DNAunzipping1.0.kamograph", "dump.DNAunzipping1.1.kamograph",
"dump.DNAunzipping1.2.kamograph", "dump.DNAunzipping1.3.kamograph", "dump.DNAunzipping1.4.kamograph",
"dump.DNAunzipping1.5.kamograph", "dump.DNAunzipping1.6.kamograph", "dump.DNAunzipping1.7.kamograph"]


averagekamos = open("average.txt", "w")
for input in inputs:
    total = 0
    frame = 0
    timestep = 0
    ndata = 0
    with open("/home/s1718990/DNAUnzipping/cos20/outputs/"+input+".txt") as fp:
        while timestep <= maxtime:
            print (timestep)
            fp.seek(0)
            for i, line in enumerate(fp):
                if (i >= (timestep/10000 * atoms/2) and i < ((timestep/10000 + 1) * atoms/2)):

                    newline = line.split()
                    if timestep >= eqtime:
                        total += int(newline[2])
                        ndata += 1
#            print(timestep,maxtime)
            timestep += 10000
            if timestep >= eqtime:
                frame += 1
        total = float(total)/float(ndata)
        averagekamos.write(str(timestep) + " " + str(total)+ " " +str(ndata) +'\n')
