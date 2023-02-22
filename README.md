# Update: New in version 0.8

Add Python Interface inside the folder Python Package
Install using PIP 
- clone this repository and go inside the FastHare folder
- `pip install ./python_package`

More details can be found in python_package/README.md

# FastHare
Source code for the paper "FastHare: Fast Hamiltonian Reduction for Large-scale Quantum Annealing, IEEE Int. Conf. on Quantum Computing and Engineering (QCE) 2022"
https://arxiv.org/abs/2205.05004

The program compresses Ising Hamiltonian and QUBOs to reduce the number of variables, and, thus, the number of Qubits for quantum codes in Quantum Annealing, Ising machine, and QAOA. The Ising Hamiltonian and QUBOs need to be first converted to Sherrington-Kirkpatrick Hamiltonians (our .net format) that do not contain any linear terms. We use the terms graph and Hamiltonian interchangably to refer to the Hamiltonians without linear terms.
 
## Compile
> make

## Usage
  
>./fasthare \<input-file in net format\> \<output-file\> $\alpha$

In each iteration, the program computes the similarity score for $n*\alpha$ edges (where $n$ is the number of nodes (variables) ).

The program will write to console in the following format

> {input_file},{n},{m},{n1},{run_time}

The program will output two files

>{output-file}_compress

The first line of this file consists of two integers $n'$ and $m'$, where $n'$ is the number of nodes and $m'$ is the number of edges in the compressed graph/Hamiltonian.

>{output-file}_flip

+ The first line contains two integers $n'$ and   $n$, where $n'$ is the number of nodes/variables *after* compression and $n$ is the number of nodes *before* compression. 
+ In the next $n'$ line, the $i$-th line contains
	+ The fisrt number: the number of nodes that are compressed to the $(i-1)$-th node. 
	+ The ids of the nodes that are compressed to the $(i-1)$-th node in the compressed graph.
+ The last line consists of $n$ number, each number is 0 or 1. If the $i$-th number is one, the corresponding value of the $(i-1)$-th qubit/variable is flipped in the graph before the compression.

## Examples

>./fasthare exp1.net exp1 1.0

>exp1.net,2,1,1,0.000114

| <span>exp1</span>.net	|exp1_compress        |exp1_flip        |
|-----------------------|---------------------|-----------------|
| 2 1 <br>0 1 3	|1 0                  |1 2 <br> 2 1 0 <br> 0 0           |
| $\displaystyle\min_{\mathbf{x}\in\{-1,+1\}^2} -[3 s_0 s_1]$	|Compressed into $n'=1$ node and $m'=0$ edges        | $s_0=s_1 = s_0'$. <br>No qubits flipped before compression. <br>The optimal solutions are either <br>$s_0=s_1 = s_0' = -1$ or <br>$s_0=s_1 = s_0' = +1$.       |



>./fasthare exp2.net exp2 1.0

>exp2.net,2,1,1,0.000254

| <span>exp2</span>.net	|exp2_compress        |exp2_flip        |
|-----------------------|---------------------|-----------------|
| 2 1 <br>0 1 -3	|1 0                  |1 2 <br> 2 1 0 <br> 0 1           |
| $\displaystyle\min_{\mathbf{x}\in\{-1,+1\}^2} -[-3 s_0 s_1]$	|Compressed into $n'=1$ node and $m'=0$ edges |  $s_0= - s_1 = s_0'$. <br>The qubit $s_0$ is not flipped and the qubit $s_1$ is *flipped* before compression. <br>The optimal solutions are either<br>$s_0=-s_1 = s_0' = -1$ or <br>$s_0= -s_1 = s_0' = +1$.  |

>./fasthare exp3.net exp3 0.2

> exp3.net,10,14,3,0.000941

Changing $\alpha=0.2$ make the program runs faster, albeit, with less effective compression

| <span>exp3</span>.net	|exp3_compress        |exp3_flip        |
|-----------------------|---------------------|-----------------|
| 10 14<br>0 1 -2<br>0 2 -4<br>1 2 2<br>1 3 -1<br>2 3 -6<br>3 4 3<br>4 5 4<br>4 6 4<br>5 6 -1<br>6 7 3<br>6 8 3<br>7 8 -1<br>7 9 3<br>8 9 2	|3 3<br>0 1 1<br>0 2 3<br>1 2 3   |3 10<br>1 8<br>2 9 7<br>7 6 5 4 3 2 0 1<br>0 1 1 0 0 0 0 0 0 0         |
| ...	|Compressed into the below Hamiltonian with $n'=3$ nodes and $m'=3$ edges.<br> $\displaystyle\min_{\mathbf{x}\in\{-1,+1\}^3} -[ s_0 s_1 + 3 s_0 s_2 + 3 s_1 s_2]$	  |  $s_8=  s_0'$<br>$s_9=  s_7=s_1'$<br>$s_6= s_5=s_4 =s_3 = - s_2 = =s_1 = s_0= s_2'$ <br>The qubits $s_1$ and $s_2$ are *flipped* before the compression.  |


## Convert MQLib instances (Max-cut) and QUBO into .net format (Hamiltonian without linear terms)
### Download and convert MQLib instances
Download the instance from AWS S3
> wget https://mqlibinstances.s3.amazonaws.com/{instance-name}.zip
> unzip {instance-name}.zip

Convert the instance into .net format
> python3 parse.py {instance-name}

### Convert from QUBO to .net file

>python3 qubo2ising.py {qubo-file} {net-file}

The program read the QUBO in {qubo-file} and write the Hamiltonian without linear terms to {net-file}
