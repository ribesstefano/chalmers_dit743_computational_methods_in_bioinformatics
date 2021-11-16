##### Author: Stefano Ribes *ribes@chalmers.se*

This PDF has been genered by compiling the `README.md` Github Markdown file.

---

# Assignment 1

The programs have been written in C99 and can be compiled by running the included Makefile:

```bash
make ; make run
```

The `global_alignment.exe` program has been modified to take in two sequences as input arguments. The default sequences are: "ATCGAT" and "ATACGT".

## Question 1+2+3: Global Alignment

By leaving the `X` and `Y` sequences to the default strings "ATCGAT" and "ATACGT" respectively,
the output alignment of `global_alignment.exe` is the following:

```
[INFO] Score matrix:
          A    T    A    C    G    T
     0   -2   -4   -6   -8  -10  -12
A   -2    2    0   -2   -4   -6   -8
T   -4    0    4    2    0   -2   -4
C   -6   -2    2    3    4    2    0
G   -8   -4    0    1    2    6    4
A  -10   -6   -2    2    0    4    5
T  -12   -8   -4    0    1    2    6

[INFO] Trace matrix:
          A    T    A    C    G    T
     0    0    0    0    0    0    0
A    0    3    2    3    2    2    2
T    0    1    3    2    2    2    3
C    0    1    1    3    3    2    2
G    0    1    1    3    3    3    2
A    0    3    1    3    3    1    3
T    0    1    3    1    3    1    3

* Alignment Sequence:

AT-CGAT
|| || |
ATACG-T

[INFO] Percent identity: 71.43%
[INFO] Hamming distance: 2
```

The percentage identity is calculated based on the total alignment length, i.e. including all the indel.

The Hamming distance is computed by getting the difference between the total alignment length and the number of matches.

## Question 4: Local Alignment

The output of the program `local_alignment.exe` is the following:

```
[INFO] Best score 6 at position (5, 9)
[INFO] Score matrix:
          H    D    A    G    A    W    G    H    E    Q
     0    0    0    0    0    0    0    0    0    0    0
P    0    0    0    0    0    0    0    0    0    0    0
A    0    0    0    2    0    2    0    0    0    0    0
W    0    0    0    0    1    0    4    2    0    0    0
H    0    2    0    0    0    0    2    3    4    2    0
E    0    0    1    0    0    0    0    1    2    6    4
A    0    0    0    3    1    2    0    0    0    4    5
E    0    0    0    1    2    0    1    0    0    2    3

[INFO] Trace matrix:
          H    D    A    G    A    W    G    H    E    Q
     0    0    0    0    0    0    0    0    0    0    0
P    0    0    0    0    0    0    0    0    0    0    0
A    0    0    0    3    0    3    0    0    0    0    0
W    0    0    0    0    3    0    3    2    0    0    0
H    0    3    0    0    0    0    1    3    3    2    0
E    0    0    3    0    0    0    0    3    3    3    2
A    0    0    0    3    2    3    0    0    0    1    3
E    0    0    0    1    3    0    3    0    0    3    3

AW-HE
|| ||
AWGHE

[INFO] Percent identity: 80.00%
[INFO] Hamming distance: 1
```


## Question 5: Levenshtein Distance

The Levenshtein distance is computed in a recursive way. The output of the program `levenshtein.exe` is as follows:

```
[INFO] Score matrix:
          A    T    A    C    G    T
     0   -2   -4   -6   -8  -10  -12
A   -2    2    0   -2   -4   -6   -8
T   -4    0    4    2    0   -2   -4
C   -6   -2    2    3    4    2    0
G   -8   -4    0    1    2    6    4
A  -10   -6   -2    2    0    4    5
T  -12   -8   -4    0    1    2    6

[INFO] Trace matrix:
          A    T    A    C    G    T
     0    0    0    0    0    0    0
A    0    3    2    3    2    2    2
T    0    1    3    2    2    2    3
C    0    1    1    3    3    2    2
G    0    1    1    3    3    3    2
A    0    3    1    3    3    1    3
T    0    1    3    1    3    1    3

AT-CGAT
|| || |
ATACG-:

[INFO] Percent identity: 71.43%
[INFO] Levenshtein dist: 2
```

Note that the implementation is not optimized.

## Question 6+7: Optimal Paths

In order to find all the optimal path, I treated the score and trace matrices together as a *directed graph*. In this way, all the paths can be found by using a backtracking-style algorithm.

The test can be run by typing:

```bash
./global_alignment.exe ATTA ATTTTA
```

The final output of the search is the following:

```
[INFO] Score matrix:
          A    T    T    T    T    A
     0   -2   -4   -6   -8  -10  -12
A   -2    2    0   -2   -4   -6   -8
T   -4    0    4    2    0   -2   -:
T   -6   -2    2    6    4    2    0
A   -8   -4    0    4    5    3    4

[INFO] Trace matrix:
          A    T    T    T    T    A
     0    0    0    0    0    0    0
A    0    3    2    2    2    2    3
T    0    1    3    3    3    3    :
T    0    1    3    3    3    3    2
A    0    3    1    1    3    3    3

* Alignment Sequence:
A--TTA
|  |||
ATTTTA
[INFO] Percent identity: 66.67:
[INFO] Hamming distance: 2
[INFO] All optimal paths:
* Path n.1: List of coordinates: (4, 6) (3, 5) (2, 4) (1, 3) (1, 2) (1, 1) (0, 0) 
* Alignment Sequence:
A--TTA
|  |||
ATTTTA
[INFO] Percent identity: 66.67:
[INFO] Hamming distance: 2

* Path n.2: List of coordinates: (4, 6) (3, 5) (2, 4) (2, 3) (1, 2) (1, 1) (0, 0) 
* Alignment Sequence:
A-T-TA
| | ||
ATTTTA
[INFO] Percent identity: 66.67:
[INFO] Hamming distance: 2

* Path n.3: List of coordinates: (4, 6) (3, 5) (2, 4) (2, 3) (2, 2) (1, 1) (0, 0) 
* Alignment Sequence:
AT--TA
||  ||
ATTTTA
[INFO] Percent identity: 66.67%
[INFO] Hamming distance: 2

* Path n.4: List of coordinates: (4, 6) (3, 5) (3, 4) (2, 3) (1, 2) (1, 1) (0, 0) 
* Alignment Sequence:
A-TT-A
| || |
ATTTTA
[INFO] Percent identity: 66.67%
[INFO] Hamming distance: 2

* Path n.5: List of coordinates: (4, 6) (3, 5) (3, 4) (2, 3) (2, 2) (1, 1) (0, 0) 
* Alignment Sequence:
AT-T-A
|| | |
ATTTTA
[INFO] Percent identity: 66.67%
[INFO] Hamming distance: 2

* Path n.6: List of coordinates: (4, 6) (3, 5) (3, 4) (3, 3) (2, 2) (1, 1) (0, 0) 
* Alignment Sequence:
ATT--A
|||  |
ATTTTA
[INFO] Percent identity: 66.67%
[INFO] Hamming distance: 2

[INFO] Number of optimal paths found: 6
```

Please refer to the functions `print_all_paths()` and `print_all_paths_util()` for the implementation of the searching algorithm.
For finding the paths, before calling `print_all_paths()`, I modified the fill matrices step by storing all the scores of all cells.

I counted the optimal paths while listing/finding them.


## Notes

## Questions

* `Minimum Score`: If I don't complete the last parts, do I still complete the assignment? Are the last parts optional?
* `Optimal Path`: What's an *optimal path*?
* `Percentage Identity`: Where can I find more information about which percentage to choose?