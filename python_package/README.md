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

>Using list of triples i, j, h_ij
>([], [0, 0], [1, -1], 1.9e-05)
>Reading triples i, j, h_ij from file
>[(0, 1, 1.0), (0, 2, 3.0), (1, 2, 3.0)]
>[2, 2, 2, 2, 2, 2, 2, 1, 0, 1]
>4.4e-05

The output is a tuple of <compressed_hamiltonian, spin_mapping, spin_sign, running_time>.
compressed_hamiltonian: triples of (i, j, h_ij)
spin_mapping: The ith entry indicates which spin in the compressed hamiltonian that the ith spin in the input Hamiltonian mapped to
spin_sign: The ith entry is -1 if the ith spin in the input Hamiltonian is reversed (+1 otherwise)
running_time: run time in second(s)


