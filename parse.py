import sys
import numpy as np
import random

k = 2
file = sys.argv[1]
out = file + ".net"
file = file + ".txt"

with open(file) as f:
    lines = f.readlines()
    
fo = open(out, "w")


n = 0
m = 0
cnt = 0
p = 0
for line in lines:
    if line.startswith("#"):
        continue
    line = line.split()
    if (len(line) < 3):
        n = int(line[0])
        m = int(line[1])        
        fo.write(str(n) + " "  + str(m) + "\n")
        continue         
    edge = [int(line[0]) - 1, int(line[1]) - 1, float(line[2])]        
    p = max(p,edge[0])
    p = max(p,edge[1])
    fo.write(str(edge[0]) + " " + str(edge[1]) + " " + str(-edge[2]) + "\n")    
    cnt = cnt + 1
fo.close()    
#print(file + "," + str (n) + "," + str(m))
 
