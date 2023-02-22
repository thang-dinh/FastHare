Python Interface for FastHare

Installation
------------
This requires Python 3.7+ and c++14

 - clone this repository
 - `pip install ./python_package`

Test
---------

```python
import fasthare as m
print("Using list of triples i, j, h_ij")
sk_ising = [(0, 1, -3)] # SK Hamiltonian: min - (-3 * x_0  * x_1 )
print(m.fasthare_reduction(sk_ising) )


# Or using the exp3.net in the "tests" folder
print("Reading triples i, j, h_ij from file")
rh, map, sign, time = m.fasthare_reduction(file ="exp3.net",alpha=0.2)
print(rh) 
print(map)
print(sign)
print(time)
```
Output

```Using list of triples i, j, h_ij
([], [0, 0], [1, -1], 1.9e-05)
Reading triples i, j, h_ij from file
[(0, 1, 1.0), (0, 2, 3.0), (1, 2, 3.0)]
[2, 2, 2, 2, 2, 2, 2, 1, 0, 1]
[1, -1, -1, 1, 1, 1, 1, 1, 1, 1]
4.4e-05
```

The output is a tuple of <compressed_hamiltonian, spin_mapping, spin_sign, running_time>.

+ compressed_hamiltonian: triples of (i, j, h_ij)
+ spin_mapping: The ith entry indicates which spin in the compressed hamiltonian that the ith spin in the input Hamiltonian mapped to
+ spin_sign: The ith entry is -1 if the ith spin in the input Hamiltonian is reversed (+1 otherwise)
+ running_time: run time in second(s)

Explanation for the output of "exp3.net"



| <span>exp3</span>.net	|compressed_hamiltonian        |spin_mapping/spin_sign        |
|-----------------------|---------------------|-----------------|
| 10 14<br>0 1 -2<br>0 2 -4<br>1 2 2<br>1 3 -1<br>2 3 -6<br>3 4 3<br>4 5 4<br>4 6 4<br>5 6 -1<br>6 7 3<br>6 8 3<br>7 8 -1<br>7 9 3<br>8 9 2	|[(0, 1, 1.0),<br> (0, 2, 3.0),<br> (1, 2, 3.0)]   |[2, 2, 2, 2, 2, 2, 2, 1, 0, 1]<br>[1, -1, -1, 1, 1, 1, 1, 1, 1, 1]         |
| ...	|Compressed into the below Hamiltonian with $n'=3$ nodes and $m'=3$ edges.<br> $\displaystyle\min_{\mathbf{x}\in\{-1,+1\}^3} -[ s_0 s_1 + 3 s_0 s_2 + 3 s_1 s_2]$	  |  $s_8=  s_0'$<br>$s_9=  s_7=s_1'$<br>$s_6= s_5=s_4 =s_3 = - s_2 = -s_1 = s_0= s_2'$ <br>The qubits $s_1$ and $s_2$ are *flipped* before the compression.  
