"""This is the package `HamGraph`, the data structure for Ising Hamiltonian
described in the paper `FastHare: Fast Hamiltonian Reduction for Large-scale Quantum Annealing`,
(https://arxiv.org/abs/2205.05004).

Copyright 2023 by Thang N. Dinh
Email: tndinh@vcu.edu
"""

import networkx as nx
import webbrowser 
import math


def log_norm(value, max_value = None):
  """Map the interval (min, max) to (0, 1) using arctan
  """
  if not max_value:
    return math.atan(value)
  return math.atan(value)/math.atan(max_value)


class   HamGraphSimple:
    # A dictionary of weighted edges {(u, v): w}
    # ,e.g., {('a', 'b'): -1.0, ('a', 'c'):2}
    graph = {}
    offset = 0
    linear_spin = '_linear_spin_'

    def __init__(self, mode = None, **kwargs):
        """Initialize a HamGraph using one of the three methods depending on the mode

        Args:
            mode (str, optional, default=None):
                `ising':    Initialized from an Ising model given by h and J parameters.
                `fasthare_ising': A list of triple (u, v, w) indicating edges (u, v) of weight w in the Hamiltonian Graph
                `networkx': Convert a networkx graph to a Hamiltonian graph 
            **kwards 
        """
        if mode == 'ising':
            self.from_ising(**kwargs)
        elif mode == 'fasthare_ising':
            self.from_fasthare_ising(**kwargs)
        elif mode == 'networkx':
            self.from_networkx(**kwargs)
        elif mode == 'netfile':
            self.get_net_file(**kwargs)
        else:
            self.offset = 0
            self.graph = {}
    
    
    def set_linear_spin(self, linear_spin):
        """Mark a different linear_spin.
        There is no checking on whether linear_spin is 
        an actual node label.
        """
        self.linear_spin = linear_spin

    def get_net_file(self, filename):      
      with open(filename,'r') as f:
        fh_ising = []
        for line in f:
          triple = [int(x) for x in line.split()]
          if len(triple) == 3: #skip the first line, containing n & m
            fh_ising.append(triple)
        self.from_fasthare_ising(fasthare_ising = fh_ising)

        
           

    def from_ising(self, h, J, offset = 0):
        """Convert Ising to HamGraph
        
        By default, the linear biases, given by h, will be 
        will be mapped to quadratic biases from 
        a special spin, called linear_spin = '_linear_spin_'
        """
        self.offset = offset
        # Note the minus. The signs in HamGraph 
        # is reversed to preserve the min-cut formulation
        # Linear biases
        self.graph =  { (self.linear_spin, key):(-value) for key, value in h.items() }   
        # Quadratic biases
        self.graph.update( { key:(-value) for key, value in J.items() } )

    
    def to_ising(self):
        """Convert HamGraph to Ising.
        By default, quadratic biases from linear_spin
        will be converted to linear biases (h)
        """
        h = {}
        J = {}
        # Note the minus. The signs in HamGraph 
        # is reversed to preserve the min-cut formulation
        for key, value in self.graph.items():
            if key[0] == self.linear_spin:
                h[ key[1] ] = -value
            elif  key[1] == self.linear_spin:
                h[ key[0] ] = -value
            else:
                J[key] = -value
        return h, J
    
    
    def get_networkx(self):
        """
        Return a networkx graph
        """
        G = nx.Graph()
        for e, w in self.graph.items():
            G.add_edge( e[0], e[1], weight=w)
        return G
    
    def from_networkx(self, graph, linear_spin = None):
        """Create a HamGraph from a networkx graph
        """
        if linear_spin is not None:
            self.linear_spin = list(graph.nodes)[0]
        else:
            self.linear_spin = linear_spin
        self.graph = {(u,v):d['weight'] for u, v, d in graph.edges(data=True)}
        offset = 0
    
    
    def to_fasthare_ising(self, pivot = None):
        """Convert HamGraph to fasthare fasthare_ising, replacing 
        variables with their 0-based indices
        """
        if pivot is None:
            pivot = self.linear_spin
        
        ul, vl = zip(*self.graph.keys())   
        var_names =  list(set(ul +vl))
        # Bring pivot to front if any
        if pivot in var_names:
            var_names.insert(0, var_names.pop(var_names.index(pivot) ) )
        n = len(var_names)
        var_map = { label:id for label, id in zip(var_names, range(n) ) }
        # Exclude any self-loops (otherwise C++ FastHare may crash)
        fh_ising =  [ (var_map[key[0]], var_map[key[1]], value) for key, value in self.graph.items() if key[0] != key[1]]
        return fh_ising, var_names
    

    def from_fasthare_ising(self, fasthare_ising, var_names = [], offset = 0, linear_spin = None):
        """Convert 0-based index fasthare_ising, a list of triples (u, v, w), to HamGraph, 
        replacing indices with variable names if provided.
        The function assumes variables' indices are within a continuous
        range starting at zero.
        """
        self.offset = offset
        if var_names is not None:
            if linear_spin is None:
              self.graph = { (i, j): hij for i, j, hij in fasthare_ising }
            else:
              self.graph = {}
              for (u, v, value) in fasthare_ising:
                if u == linear_spin:
                  self.graph[(self.linear_spin, v)] = value
                elif v == linear_spin:
                  self.graph[(self.linear_spin, u)] = value
                else:
                  self.graph[(u, v)] = value
        else:
            if linear_spin is None:
              self.graph = { (var_names[i], var_names[j]): hij for i, j, hij in fasthare_ising }
            else:
              self.graph = {}
              for (u, v, value) in fasthare_ising:
                if u == linear_spin:
                  self.graph[(self.linear_spin, var_names[v])] = value
                elif v == linear_spin:
                  self.graph[(self.linear_spin, var_names[u])] = value
                else:
                  self.graph[(var_names[u], var_names[v])] = value
        