import neal
from FastHareComposite import FastHareComposite

# A toy Ising Hamiltonian
h= {0: -3, 1: -4, 2: 4, 3: 4} 
J= {(0, 2): 1, (0, 3): -4, (1, 3): -2, (2, 3): 3}
# Use FastHare composite to preprocess instance before sovling
fh_sa = FastHareComposite(sa)
# Solve the compressed Ising using simulated annealing (SA)
sample_set_fh = fh_sa.sample_ising(h, J, num_reads = 100)
# Print out the samples
print(sample_set_fh.aggregate())
print("Best solution has a minimum energy ", sample_set_fh.first.energy)