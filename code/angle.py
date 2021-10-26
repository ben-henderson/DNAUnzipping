import math
import os

lines =[]
results = []
timestep = 0
atoms = 400
maxtime = 4001000
input = "dump.DNAunzipping1.0"
with open("/home/ben/github/DNAUnzipping/outputs/"+input) as fp:
    while timestep <= maxtime:
        for i, line in enumerate(fp):
            if i >= (timestep*(atoms+9))+9 and i < (timestep+1)*(atoms+9):
                newline = line.split()
                lines.append(newline)
            else:
                break
            #lines.sort(key=lamda x: x[1])
        for x in range(int((atoms/2)-2)):
            dy = lines[(x*2)+1][3] - lines[x*2][3]
            dx = lines[(x*2)+1][2] - lines[x*2][2]
            angle = math.atan2(dy/dx)
            if angle < 2.5:
                result = [timestep*1000, lines[x*2][1], 0]

            else:
                result = [timestep*1000, lines[x*2][1], 1]

            results.append(result)
            timestep =+ 1000


print (results)
