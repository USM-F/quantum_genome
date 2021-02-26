## Genome assembly tutorial
Code for the paper [Genome assembly using quantum and quantum-inspired annealing](https://arxiv.org/abs/2004.06719)

### Usage
Example for syntetic data

```python
import numpy as np
import networkx as nx
import re

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave.system import LeapHybridSampler

import DeBruijnDNA
import AcyclicGraphDNA

API_KEY = "D-WAVE_OCEAN_SDK_API_KEY"

sequence = 'CATACACCTAA'
kmer_len, suffix_len = 3, 2

adjacency_matrix, node_labels = DeBruijnDNA.make_debr(sequence, kmer_len=kmer_len, suffix_len=suffix_len)
Q = DeBruijnDNA.to_qubo(adjacency_matrix)

q = {}
size = len(Q)
for i in range(size):
    for j in range(size):
        q[(i,j)] = Q[i][j]
            
sampler = LeapHybridSampler(token=API_KEY)
response = sampler.sample_qubo(q)

result = []
for sample, energy in response.data(['sample', 'energy']): 
    result.append(sample)
    result.append(energy)

result.append(response.info) # view timings
spins, energy = [result[0][i] for i in result[0].keys()], result[1]
```
In spins equal 1 is coded node position i = (spin index) % (sequence length - 2)

For more detailed examples see example.ipynb
