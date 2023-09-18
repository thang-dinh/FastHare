
"""This is the package `FastHareComposite`.
It provides a simple composite to work with the `Samplers` from  dwave-ocean-sdk.

Copyright 2023 by Thang N. Dinh
Email: tndinh@vcu.edu
"""

from dimod.core.composite import Composite
from dimod.core.sampler import Sampler
from dimod.sampleset import SampleSet, concatenate
from dimod.vartypes import Vartype
import math
from HamGraph import HamGraph
import dimod
import fasthare as fh

__all__ = ['FastHareComposite']

class FastHareComposite(Sampler, Composite):
    """FastHare reduction of Hamiltonian.
    """
    children = None
    parameters = None
    properties = None

    def __init__(self, child):
        self.children = [child]

        self.parameters = parameters = {'alpha': None}
        parameters.update(child.parameters)

        self.properties = {'child_properties': child.properties}

    def sample(self, bqm, *, alpha=None, verbose=True, **kwargs):
        """Sample from the binary quadratic model.

        Args:
            bqm (:class:`~dimod.BinaryQuadraticModel`):
                Binary quadratic model to be sampled from.

            alpha (float, optional, default=1.0):
                How much effort in trying to compress (time complexity O(\alpha n^2)).
            verbose (bool, optional, default=True):
                Print out reduction efficiency

        Returns:
            :class:`.SampleSet`
        
        Examples:
            This example runs FastHare applied to one variable of an Ising problem.

        """
        #assert(bqm.vartype in [dimod.SPIN, dimod.BINARY], "Unsupported VarType for variables. FastHare only supports SPIN and BINARY")
        fh_bqm = bqm.copy()
        # print("Bqm model: ", bqm)
        h, J, offset = fh_bqm.to_ising() 
        #print("Original Ising: ", h, J)
        #Remove zero coefficients/biases
        hp = {key:value for key, value in h.items() if abs(value) > 1e-5}
        #Remove self-loops
        Jp = {key:value for key, value in J.items() if abs(value) > 1e-5 and key[0] != key[1]}
        # print(hp, Jp)
        hg = HamGraph("ising", h = hp, J = Jp, offset = offset)
        fh_ising, var_names = hg.to_fasthare_ising()
        # Make a set of ``ghost'' variables that got eliminated
        sv = set(var_names)
        ghost_var = [v for v in bqm.variables if not v in sv]
        # print("Stranded variables ", ghost_var)
        # Set alpha to log number of variables
        if not alpha:
            self.alpha = int(math.log(len(bqm.variables)))+1
        else:
            self.alpha = alpha
        # print("Input Ising for FastHare: ", fh_ising)
        fh_reduced_ising, fh_map, fh_sign, runtime = fh.fasthare_reduction(sk_ising = fh_ising, alpha = self.alpha )
        # print("reduced map/sign \n",fh_map,"\n", fh_sign)
        fbqm = lambda x: x
        if bqm.vartype == dimod.BINARY:
            fbqm = lambda x: (x+1)/2.0
        if not fh_reduced_ising:           
           if verbose:
             print("FastHare: Reduced 100% variables. Exact solution found.")
           sample =  {var_names[i]:fbqm(fh_sign[i]*fh_sign[0]) for i in range(1, len(var_names) )} 
           sample.update({v:1 for v in ghost_var})
        #    print("Optimal solution found ", sample)
           return SampleSet.from_samples_bqm(sample, bqm)
        else:
            hgr = HamGraph("fasthare_ising", fasthare_ising = fh_reduced_ising, linear_spin = 0)
            hr, Jr = hgr.to_ising()
            
            bqm_r = dimod.BinaryQuadraticModel(hr, Jr, dimod.SPIN)
            if verbose:
                print("FastHare: Reduced {}/{} variables ({:.2f}% reduction).".format(\
                    len(bqm.variables) - len(bqm_r.variables),\
                    len(bqm.variables),\
                    100.0*(len(bqm.variables) - len(bqm_r.variables))/ len(bqm.variables)))
            sampleset_r = self.child.sample(bqm_r, **kwargs)
            # print(sampleset_r)
            # If there is no linear_spin, then include 0th variable in the solution
            b_index = 1
            if (var_names[0] != hg.linear_spin):
              b_index = 0
            sampleset = [] 
            for sample in sampleset_r:
                f0 = lambda idx: 1 if idx == 0 else sample[idx]  
                sign_h = f0(fh_map[0]) * fh_sign[0]
                sol = { var_names[i]: fbqm( f0( fh_map[i] ) * fh_sign[i] * sign_h) for i in range(b_index,  len(var_names))}    
                #print(sol)
                sol.update({v:1 for v in ghost_var})
                sampleset.append( sol )
            return SampleSet.from_samples_bqm(sampleset, bqm)
                

    
