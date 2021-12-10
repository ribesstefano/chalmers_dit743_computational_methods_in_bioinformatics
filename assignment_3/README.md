# Author: Stefano Ribes

# Assignment 3: Main chain tracing

This PDF document can be generated from the markdown file `README.md`.

## Question 1

### Setup

The first part of the assignment has been written in Python (version 3) and can be found in file `main_chain.py`. The script requires the `matplotlib` package in order to visualize the resulting chain.

The program can be run with the following command:
```bash
python3 main_chain.py
```

### Algorithm

The algorithm starts by generating all the possible atom pairs wich have a distance below a certain threshold.
I started with a threshold of 3.8 Å, but the found pairs didn't cover all atoms in the file.
Because of that, I set the threshold to 4.5 Å instead. In this way, I'm able to include all atoms available in the pairs.

Once all the proper pairs of atoms have been found, the algorithm starts building the alpha Carbon chain. The first step consists of setting the pairs as a pool of possible chain segments.
A first pair/segment is extracted from the pool. From that, the algorithm will compare the extremes of the chain with the remaining pairs in the pool.

In particular, if an atom in a pair is matching an extreme, the pair is appended to the chain. An example follows:
```
# Chain starts from atom pair: (2, 4)

chain: [2, 4] + atom (10, 9) -> chain: [2, 4]         # Neither atom 10 nor 9 are added
chain: [2, 4] + atom (4, 1) -> chain: [2, 4, 1]       # Atom 1 is added to the chain
chain: [2, 4, 1] + atom (3, 1) -> chain: [2, 4, 1, 3] # Atom 3 is added to the chain
...
```

The final output of the program is:

```
[DEBUG] Appending (3, 5), chain: [2, 4]
[DEBUG] Appending (3, 1), chain: [2, 4]
[DEBUG] Appending (10, 8), chain: [2, 4]
[DEBUG] Appending (10, 9), chain: [2, 4]
[DEBUG] Appending (4, 1), chain: [2, 4, 1]
[DEBUG] Appending (7, 6), chain: [2, 4, 1]
[DEBUG] Appending (7, 9), chain: [2, 4, 1]
[DEBUG] Appending (3, 5), chain: [2, 4, 1]
[DEBUG] Appending (3, 1), chain: [2, 4, 1, 3]
[DEBUG] Appending (10, 9), chain: [2, 4, 1, 3]
[DEBUG] Appending (5, 6), chain: [2, 4, 1, 3]
[DEBUG] Appending (7, 6), chain: [2, 4, 1, 3]
[DEBUG] Appending (7, 9), chain: [2, 4, 1, 3]
[DEBUG] Appending (3, 5), chain: [2, 4, 1, 3, 5]
[DEBUG] Appending (10, 9), chain: [2, 4, 1, 3, 5]
[DEBUG] Appending (5, 6), chain: [2, 4, 1, 3, 5, 6]
[DEBUG] Appending (7, 9), chain: [2, 4, 1, 3, 5, 6]
[DEBUG] Appending (10, 8), chain: [2, 4, 1, 3, 5, 6]
[DEBUG] Appending (10, 9), chain: [2, 4, 1, 3, 5, 6]
[DEBUG] Appending (7, 6), chain: [2, 4, 1, 3, 5, 6, 7]
[DEBUG] Appending (10, 8), chain: [2, 4, 1, 3, 5, 6, 7]
[DEBUG] Appending (10, 9), chain: [2, 4, 1, 3, 5, 6, 7]
[DEBUG] Appending (7, 9), chain: [2, 4, 1, 3, 5, 6, 7, 9]
[DEBUG] Appending (10, 8), chain: [2, 4, 1, 3, 5, 6, 7, 9]
[DEBUG] Appending (10, 9), chain: [2, 4, 1, 3, 5, 6, 7, 9, 10]
[DEBUG] Appending (10, 8), chain: [2, 4, 1, 3, 5, 6, 7, 9, 10, 8]
chain: 8 10 9 7 6 5 3 1 4 2
gold:  8 10 9 7 6 5 3 1 4 2
```

## Question 2

### Setup

The second part of the assignment has been written in Python (version 3) and can be found in file `main_chain_free_atoms.py`.

The program can be run with the following command:
```bash
python3 main_chain_free_atoms.py
```

### Algorithm

The algorithm isn't complete. A function for calculating 3D angles has been written.