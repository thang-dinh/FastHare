import dimod
import sys
from collections import defaultdict
Q = defaultdict(float)
file = sys.argv[1]
with open(file) as f:
    lines = f.readlines()
out = sys.argv[2]
fo = open(out, "w")
for line in lines:
    if line.startswith("#"):
        continue
    line = line.split()
    if (len(line) < 3):
        n = int(line[0])
        m = int(line[1])                
        continue         
    edge = [int(line[0]), int(line[1]), float(line[2])]        
    u = edge[0]
    v = edge[1]
    w = edge[2]    
    Q[(u,v)] += w
#print(Q)
K = dimod.qubo_to_ising(Q, 0)
n = n + 1
m = len(K[0]) + len (K[1])
fo.write(str(n)+" " +str(m)+"\n")
for i in K[0]:
    fo.write("0 " + str(i) + " " + str(K[0][i])+"\n")
for (u,v) in K[1]:
    fo.write(str(u) + " " + str(v) + " " + str(K[1][(u,v)])+"\n")