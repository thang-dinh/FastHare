"""This is the package `HamGraph`, the data structure for Ising Hamiltonian
described in the paper `FastHare: Fast Hamiltonian Reduction for Large-scale Quantum Annealing`,
(https://arxiv.org/abs/2205.05004).

Copyright 2023 by Thang N. Dinh
Email: tndinh@vcu.edu
"""

import networkx as nx
import webbrowser 
import math
import os
import HamGraph as hg
import IPython
import plotly.express as px
import gravis as gv



def log_norm(value, max_value = None):
  """Map the interval (min, max) to (0, 1) using arctan
  """
  if not max_value:
    return math.atan(value)
  return math.atan(value)/math.atan(max_value)

def draw_gv(G, html_file='netvis.html', display_method = None, sample = None):
  """Visualize a networkx graph using Gravis library

  Args:
    html_file (str, optional, default = 'netvis.html'):
        Name of the html file to write the graph to
    display_method (str, optional, default = None):
        One of the three values 'notebook', 'browser'. Default is None
        Setting display_method to 'notebook' will attempt to display the 
        html file using IPython within Jupyter notebook and such. The 'browser' 
        will attempt to open the html file in the browser
    sample (dictionary, optional, default = None):
        A dictionary contains +1/-1 classification of nodes
  """
  color_map = 'RdBu'
  # Compute pagerank after ignoring the signs on edge weights
  G2 = G.copy()
  for u, v, a in G2.edges(data=True):
     a['weight'] = abs(a['weight'])
  # Set nodes' color
  if sample:
     for u in G.nodes():
        if u in sample:
            if sample[u] == -1:
                G.nodes[u]['color'] = "#3399CC" #
            elif sample[u] == 1:
                G.nodes[u]['color'] = "#CDB99C" #
  # Set hover to be a list of neighbors + weight
  for u in G.nodes():
     nbw = { v: G[u][v]['weight'] for v in G[u]}
     sorted_nb = list( sorted(nbw.items(), key = lambda item: -abs(item[1])))
     G.nodes[u]['hover'] = "Node:"+ str(u)+ ", deg: "+str(len(nbw)) +"<br/>"+\
        "<br>".join( [ "{:>5}| {:>10}".format(v, w) for v, w in  sorted_nb])
            
  # print(nx.get_node_attributes(G, "color"))
  # Set the size to be the pagerank
  pr = nx.pagerank(G2, alpha = 0.9)
  nx.set_node_attributes(G, pr, 'size')

  # Edge width = abs( edge weight)
  for u, v, a in G.edges(data=True):
     a['width'] = abs(a['weight'])
     if u == 'h' or v == 'h':
        a['opacity'] = 0.1 
  # Set edge colors
  ew = nx.get_edge_attributes(G, "weight")
  minw, maxw = min(ew.values()), max(ew.values())
  for u,v,a in G.edges(data=True):
    if ('weight' in  a):
      cl = 0
      if a['weight'] >= 0:
        cl = log_norm(a['weight'], maxw)
      else:
        cl = log_norm(a['weight'], -minw)
      color = px.colors.sample_colorscale(color_map, [cl*0.4 + 0.5])
      a['color'] = color[0]
      a['hover'] = str(a['weight'])
                # + " - " + str(a['color']) + " - "+str(a['value'])

  # Display   
  fig = gv.d3(G, use_node_size_normalization=True, node_size_normalization_max=30,
      use_edge_size_normalization=True, edge_size_data_source='width',show_menu=True,\
        edge_curvature=0.25, edge_hover_tooltip=True, graph_height=700)
  if os.path.isfile(html_file):
    print('Warning: File exists. Overwritten ', html_file)
    os.remove(html_file)
  fig.export_html(html_file)
  if display_method == 'notebook':
    IPython.display.display(IPython.display.HTML(filename=html_file))
  elif display_method == 'browser':
    path = 'file:///'+os.getcwd()+'/'+html_file
    print('Opening in browser...', path)
    hnd = webbrowser.open(path)

class   HamGraph:
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
        else:
            self.offset = 0
            self.graph = {}
    
    
    def set_linear_spin(self, linear_spin):
        """Mark a different linear_spin.
        There is no checking on whether linear_spin is 
        an actual node label.
        """
        self.linear_spin = linear_spin

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
        if not linear_spin:
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
        if not var_names:
            if not linear_spin:
                self.graph = { (i, j): hij for i, j, hij in fasthare_ising }
            else:
                self.graph = {}
                for key, value in fasthare_ising.items():
                   if key[0] == linear_spin:
                      self.graph[(self.linear_spin, key[1])] = value
                   elif key[1] == linear_spin:
                      self.graph[(self.linear_spin, key[0])] = value
                   else:
                      self.graph[key] = value
        else:
            if not linear_spin:
                self.graph = {}
                for key, value in fasthare_ising.items():
                   if key[0] == linear_spin:
                      self.graph[(self.linear_spin, var_names[key[1]])] = value
                   elif key[1] == linear_spin:
                      self.graph[(self.linear_spin, var_names[key[0]])] = value
                   else:
                      self.graph[(var_names[key[0]], var_names[key[1]])] = value
        
    
    def draw(self, html_file='netvis.html', display_method = None, sample = None):
       """Visualize the Hamiltoninan Graph using Gravis libary

       Args:
            html_file (str, optional, default = 'netvis.html'):
                The name of the html file. Default value: `netvis.html'
            display_method (str, optional, default = None):
                One of the three values 'notebook', 'browser'. Default is None
                Setting display_method to 'notebook' will attempt to display the 
                html file using IPython within Jupyter notebook and such. The 'browser' 
                will attempt to open the html file in the browser
            sample (dictionary, optional, default = None):
                A dictionary contains +1/-1 classification of nodes
       """
       G = self.get_networkx()
       # Replace the linear_spin with "h"
       G = nx.relabel_nodes(G, {self.linear_spin:"h"}, copy=False)
       if "h" in G.nodes():
          G.nodes["h"]['color'] = '#ffffff'
          G.nodes["h"]["border_color"]='#333333'
          G.nodes["h"]["border_size"]=2

    
       
       draw_gv(G, html_file, display_method, sample)

       

            
