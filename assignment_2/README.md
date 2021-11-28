# Assignment 2: Domain assignment

## Purpose

The C programs `pdb_io.c`, `atom_array.c` and `residue_array.c`  demonstrate file handling and input/output.  

## Description

The program `pdb_io` takes the name of a Protein Data Bank (PDB) file as a
command line argument, reads each line of that file, splits up each ATOM line
into individual fields (placing their values in separate variables), then
writes out the contents of each ATOM line in three ways:

	1. (most of) the individual fields, one field per line;
	2. the single string containing the whole ATOM line (as read);
	3. (most of) the individual fields, in PDB format.

The program `atom_array` first reads PDB atom records into an array of
"atom" structures, then writes out the content of this array in PDB format.

The program `residue_array` first reads PDB atom records into an array of
"residue" structures, then writes out the content of this array in PDB format.

Later in this course you will write programs that have to read/write
Protein Data Bank files, and you might find it useful to use these
programs as a starting point.


These programs can be compiled using the Makefile (see below).
You don't *need* to understand (or even look at) the contents of the
Makefile, or how the `make` utility works.  But, if you are interested,
there are some explanatory comments in the Makefile.


## Instructions

Copy the files from this directory into a directory of your own.
I recommend that you use a new directory for each practical.

To compile the programs, type:

```bash
make
```

To run the programs with data file '1CRN.pdb', type:

```bash
./pdb_io 1CRN.pdb
./atom_array 1CRN.pdb
./residue_array 1CRN.pdb
```

To remove those files that can be recompiled from the source code, type:

```bash
make clean
```

## Major modifications to the original C code

List of modifications:

* Cleaned up formatting and style
* Created separate header files with all data structure definitions and methods
* Removed all global variables
* In `atom_array.c`, array indexing now starts from zero
* Written generic `read_data` function: it now accepts a function pointer to a callback function. Such callback expects the PDB line index, the read PDB entry and any user data. The prototype of the callback function looks something like the following:

```c
typedef void (*callback_ptr)(const PdbEntry*, int*, void* user_data); // Type definition

void atom_callback(const PdbEntry* entry, int* line_idx, void* user_data); // Example of definition

int main(int argc, char **argv) {
  int numAtoms, i;
  if (argc < 2) {
    (void) fprintf(stderr, "usage: atom_array file.pdb\n");
    exit(0);
  }
  numAtoms = read_data(argv[1], &atom_callback, (void*)atom);
  for (i = 0; i < numAtoms; ++i) {
    write_pdb_atom(
      atom[i].serial,
      atom[i].atomName,
      atom[i].altLoc,
      atom[i].resName,
      atom[i].chainID,
      atom[i].resSeq,
      atom[i].iCode,
      atom[i].centre);
  }
  return 0;
}
```
## Questions and Outputs

### Question 1: Distance Map Generation

The code of this question can be found in file `make_distance_map.c`.


### Question 2: DOMAK Partition

The code of this question can be found in file `domak_partition.c`.
The code implements the naïve algorithm for two segments partioning described in *Continuous and discontinuous domains: An algorithm for the automatic generation of reliable protein domain definitions*, by S. Siddiqui and J. Barton, 1995.

The implementation only stores *heavy atoms*, *i.e.* alpha carbon atoms, into an array for the comparisons. In case a residue doesn't contain any of such atoms, the program will terminate and output an error message. This means that we are assuming there's at least one atom per residue to compare with the others.

In order to speedup computation, the atom distances are stored into a lookup table once computed.

Please note that the implementation is not optimized for space, but rather for execution time.

#### Output of Running on 1CDH
```
[INFO] Maximum split value: 26804.250, corresponding to index 97.
[INFO] Bar plot with normalized values:
    2 ** 2.44
    3 * 1.88
    4 * 1.50
    5 * 1.24
    6 * 1.05
    7 * 1.24
    8 * 1.05
    9  0.71
   10  0.79
   11  0.71
   12  0.77
   13  0.69
   14  0.74
   15  0.66
   16  0.71
   17  0.63
   18  0.67
   19  0.71
   20  0.74
   21  0.77
   22  0.60
   23  0.43
   24  0.45
   25  0.46
   26  0.43
   27  0.40
   28  0.37
   29  0.35
   30  0.33
   31  0.34
   32  0.35
   33  0.36
   34  0.36
   35  0.42
   36  0.42
   37  0.40
   38  0.45
   39  0.42
   40  0.43
   41  0.43
   42  0.44
   43  0.45
   44  0.51
   45  0.51
   46  0.58
   47  0.59
   48  0.60
   49  0.55
   50  0.56
   51  0.56
   52  0.63
   53  0.64
   54  0.65
   55  0.59
   56  0.55
   57  0.50
   58  0.51
   59  0.51
   60  0.52
   61  0.52
   62  0.49
   63  0.59
   64  0.60
   65  0.82
   66 * 1.03
   67 * 1.04
   68 * 1.18
   69 * 1.35
   70 * 1.81
   71 ** 2.12
   72 ** 2.52
   73 *** 3.03
   74 *** 3.71
   75 *** 3.73
   76 *** 3.07
   77 *** 3.08
   78 *** 3.09
   79 ** 2.20
   80 ** 2.20
   81 * 1.89
   82 ** 2.21
   83 * 1.90
   84 * 1.90
   85 * 1.65
   86 * 1.91
   87 * 1.91
   88 * 1.91
   89 * 1.91
   90 ** 2.23
   91 ** 2.64
   92 *** 3.16
   93 *** 3.84
   94 ****** 6.06
   95 ********** 10.90
   96 ************************ 24.76
   97 **************************************************************************************************** 100.00
   98 *************************************************************************************************** 99.81
   99 *************************************************************************************************** 99.46
  100 ******************************************* 43.79
  101 ************************ 24.40
  102 *************** 15.46
  103 *************** 15.41
  104 *************** 15.36
  105 *************** 15.30
  106 ********** 10.51
  107 ******* 7.63
  108 ******* 7.60
  109 ***** 5.74
  110 **** 4.48
  111 ** 2.93
  112 ** 2.92
  113 ** 2.42
  114 ** 2.40
  115 ** 2.02
  116 ** 2.82
  117 ** 2.33
  118 ** 2.31
  119 ** 2.74
  120 ** 2.71
  121 ** 2.69
  122 ** 2.67
  123 ** 2.64
  124 ** 2.62
  125 ** 2.60
  126 ** 2.57
  127 ** 2.11
  128 * 1.76
  129 * 1.48
  130 * 1.26
  131 * 1.08
  132  0.93
  133  0.92
  134 * 1.03
  135 * 1.01
  136  1.00
  137 * 1.12
  138 * 1.11
  139 * 1.25
  140 * 1.23
  141 * 1.21
  142 * 1.38
  143 * 1.60
  144 * 1.57
  145 * 1.84
  146 * 1.80
  147 ** 2.14
  148 ** 2.59
  149 **** 4.21
  150 *** 3.08
  151 ** 3.00
  152 *** 3.82
  153 *** 3.72
  154 ** 2.70
  155 ** 2.61
  156 * 1.94
  157 * 1.48
  158 * 1.14
  159 * 1.06
  160  0.82
  161  0.93
  162  0.84
  163  0.79
  164  0.74
  165  0.69
  166  0.64
  167  0.71
  168  0.64
  169  0.71
  170  0.80
  171  0.91
  172 * 1.05
  173 * 1.24
  174 * 1.50
  175 * 1.88
  176 ** 2.44
  177 ** 2.46
```

#### Output of Running on 2CSN
```
[INFO] Maximum split value: 191464.000, corresponding to index 210.
[INFO] Bar plot with normalized values:
    2  0.51
    3  0.90
    4 * 1.28
    5 * 1.66
    6 ** 2.03
    7 * 1.07
    8 *** 3.02
    9 * 1.50
   10  0.94
   11 * 1.03
   12  0.71
   13  0.77
   14  0.57
   15  0.45
   16  0.48
   17  0.50
   18  0.53
   19  0.43
   20  0.61
   21  0.65
   22  0.95
   23  0.72
   24 * 1.54
   25 * 1.10
   26 * 1.16
   27 * 1.20
   28 * 1.23
   29 * 1.27
   30 * 1.91
   31 * 1.96
   32 *** 3.19
   33 * 1.44
   34 * 1.50
   35 * 1.12
   36 * 1.58
   37 * 1.18
   38 * 1.67
   39 * 1.24
   40 * 1.27
   41 * 1.29
   42 * 1.31
   43 * 1.33
   44 * 1.35
   45 * 1.38
   46 * 1.40
   47 * 1.42
   48 * 1.44
   49 * 1.46
   50 * 1.48
   51 * 1.50
   52 * 1.52
   53 * 1.18
   54 * 1.19
   55  0.95
   56 * 1.23
   57 * 1.24
   58 * 1.66
   59 * 1.68
   60 * 1.70
   61 * 1.72
   62 * 1.32
   63 * 1.34
   64 * 1.35
   65 * 1.36
   66 * 1.08
   67 * 1.09
   68 * 1.10
   69  0.90
   70  0.91
   71  0.76
   72  0.76
   73  0.77
   74  0.77
   75  0.95
   76 * 1.19
   77 * 1.53
   78 ** 2.03
   79 ** 2.80
   80 **** 4.09
   81 ****** 6.49
   82 ************************** 26.47
   83 ************************** 26.63
   84 **** 4.25
   85 **** 4.27
   86 ** 2.98
   87 ** 2.99
   88 *** 3.01
   89 *** 3.02
   90 ** 2.23
   91 * 1.71
   92 * 1.72
   93 ** 2.27
   94 ** 2.28
   95 ** 2.29
   96 *** 3.14
   97 ** 2.31
   98 ** 2.32
   99 ** 2.33
  100 ** 2.34
  101 ** 2.35
  102 ** 2.36
  103 ** 2.37
  104 * 1.81
  105 * 1.82
  106 * 1.82
  107 * 1.83
  108 * 1.84
  109 * 1.85
  110 ** 2.43
  111 ** 2.44
  112 ** 2.44
  113 ** 2.45
  114 ** 2.46
  115 * 1.88
  116 * 1.88
  117 * 1.89
  118 * 1.49
  119 * 1.49
  120 * 1.50
  121 * 1.91
  122 * 1.91
  123 * 1.92
  124 * 1.51
  125 * 1.51
  126 * 1.52
  127 * 1.52
  128 * 1.52
  129 * 1.52
  130 * 1.52
  131 * 1.52
  132 * 1.23
  133 * 1.53
  134 * 1.53
  135 * 1.53
  136 * 1.95
  137 * 1.95
  138 * 1.95
  139 * 1.95
  140 * 1.95
  141 * 1.95
  142 ** 2.55
  143 ** 2.55
  144 ** 2.55
  145 *** 3.49
  146 ***** 5.04
  147 ******* 7.90
  148 ******* 7.90
  149 ******* 7.89
  150 ******* 7.89
  151 ******* 7.88
  152 ******* 7.87
  153 ************** 14.03
  154 ************** 14.02
  155 ******************************* 31.62
  156 ******************************* 31.58
  157 ************* 13.95
  158 ************* 13.93
  159 ************* 13.91
  160 ************* 13.89
  161 ************* 13.87
  162 ************* 13.84
  163 ************* 13.82
  164 ******************************* 31.14
  165 ******************************* 31.08
  166 ******************************* 31.02
  167 ****************************** 30.95
  168 ****************************** 30.88
  169 ************* 13.63
  170 ************* 13.59
  171 ******* 7.59
  172 ******* 7.57
  173 ******* 7.55
  174 ******* 7.53
  175 ******* 7.50
  176 ******* 7.48
  177 ******* 7.46
  178 ******* 7.44
  179 ******* 7.41
  180 ******* 7.39
  181 ******* 7.36
  182 ******* 7.34
  183 ******* 7.31
  184 ******* 7.28
  185 ******* 7.25
  186 ******* 7.23
  187 ************ 12.84
  188 **************************** 28.85
  189 **************************** 28.73
  190 **************************** 28.61
  191 **************************** 28.48
  192 **************************** 28.35
  193 **************************** 28.22
  194 **************************** 28.09
  195 *************************** 27.95
  196 *************************** 27.82
  197 ************ 12.22
  198 ************ 12.16
  199 ****** 6.76
  200 ****** 6.70
  201 ****** 6.66
  202 ****** 6.59
  203 *********** 11.68
  204 ****** 6.49
  205 ****** 6.42
  206 ****** 6.37
  207 *********** 11.29
  208 *********** 11.21
  209 ************************* 25.11
  210 **************************************************************************************************** 100.00
  211 *************************************************************************************************** 99.27
  212 ************************************************************************************************** 98.52
  213 ************************************************************************************************* 97.77
  214 ************************************************************************************************* 97.01
  215 ************************************************************************************************ 96.24
  216 *********************************************************************************************** 95.46
  217 ********************************************************************************************** 94.67
  218 ********************************************************************************************* 93.87
  219 ********************************************************************************************* 93.06
  220 ******************************************************************************************** 92.24
  221 ******************************************************************************************* 91.42
  222 ****************************************************************************************** 90.58
  223 ********************** 22.23
  224 ********************** 22.02
  225 ********************* 21.80
  226 ************************************************************************************** 86.55
  227 ************************************************************************************* 85.66
  228 ************************************************************************************ 84.76
  229 *********************************************************************************** 83.85
  230 ********************************************************************************** 82.93
  231 ******************** 20.29
  232 ******************** 20.06
  233 ******************* 19.82
  234 ****************************************************************************** 78.53
  235 ***************************************************************************** 77.55
  236 **************************************************************************** 76.57
  237 ****************** 18.68
  238 ****************** 18.43
  239 ******* 7.98
  240 ***************** 17.75
  241 ***************** 17.49
  242 ********************************************************************* 69.09
  243 ******************************************************************** 68.03
  244 ****************************************************************** 66.97
  245 ***************************************************************** 65.89
  246 **************************************************************** 64.80
  247 *************** 15.70
  248 *************** 15.43
  249 *************** 15.15
  250 *********************************************************** 59.60
  251 ********************************************************** 58.46
  252 ********************************************************* 57.31
  253 ******************************************************** 56.15
  254 ****************************************************** 54.98
  255 ***************************************************** 53.80
  256 **************************************************** 52.62
  257 *************************************************** 51.42
  258 ************************************************** 50.21
  259 ************************************************ 49.00
  260 *********************************************** 47.77
  261 ********************************************** 46.54
  262 ********************************************* 45.29
  263 ******************************************** 44.04
  264 ****************************************** 42.78
  265 ***************************************** 41.50
  266 **************************************** 40.22
  267 ********* 9.50
  268 ********* 9.17
  269 ******** 8.84
  270 ********************************** 34.12
  271 ******* 7.95
  272 ******* 7.62
  273 ******* 7.28
  274 *************************** 27.81
  275 ************************** 26.43
  276 ************************* 25.04
  277 *********************** 23.64
  278 ***** 5.31
  279 **** 4.96
  280 **** 4.60
  281 **************** 16.99
  282 *************** 15.54
  283 ************** 14.08
  284 ************ 12.61
  285 *********** 11.13
  286 ********* 9.65
  287 ******** 8.15
  288 ****** 6.64
  289 ***** 5.12
  290 *** 3.60
  291 ** 2.06
  292  0.52
```

### Question 4: Multi Domain Partition

The code of this question can be found in file `multi_domak_partition.c`.
The code implements the methods for multi segments partioning described in *Continuous and discontinuous domains: An algorithm for the automatic generation of reliable protein domain definitions*, by S. Siddiqui and J. Barton, 1995.

In particular, the algorithm will first start assuming a single domain, which is shared list (actually a static C array). Then call the `single_segment_scan()` function, which will update the domain list with none, one or two more domains.
Afterwards, the algorithm will repeat the scan over the found domains, eventually calling `two_segment_scan_of_two_segment_domain()` when encountering a two-segments domain.
Finally, when all the domains' segments are shorter then a predefined threshold or when no new domains are found, the algorithm terminates.

#### Implementation Details

Some of the optimizations described in the paper are not implemented. The ones implemented are the following:
  * Only considering *heavy atoms*, *i.e.* alpha carbon atoms, in the residues
  * All domains assume a maximum of two segments
  * For single-segment domains, only Method 1 (*single-segment scan*) is applied
  * Minimum segment length: segments below Minimum Domain Size (MDS) are not considered valid and so not splitted (other constraints like MNCC or MSS are not considered)

A lot of segment comparisons are needed and having a lookup table of the distances is not enough for computing all the interior and exterior contacts.
A *multithreaded approach* can speed up the computation, since many segments can be computed independently. However, if the lookup table of the distances is kept shared, it would create a bottleneck among the threads since all threads would waste time requesting the lock to update it.

Because of the above issue, my implementation for partitioning exploits OpenMP and avoids using a lookup table.

#### Outputs

The protein described in PDB file `4GAF_B.pdb` can be divided in the following 7 domains:
```
[INFO] Found 7 domains:
Domain n.1
        segment n.1: (1, 42)
Domain n.2
        segment n.1: (42, 83)
Domain n.3
        segment n.1: (83, 124)
Domain n.4
        segment n.1: (124, 165)
Domain n.5
        segment n.1: (165, 206)
Domain n.6
        segment n.1: (206, 247)
Domain n.7
        segment n.1: (247, 306)
```

The protein described in PDB file `1HZH_H.pdb` can be divided in the following 9 domains:
```
[INFO] Found 9 domains:
Domain n.1
        segment n.1: (1, 59)
Domain n.2
        segment n.1: (59, 100)
Domain n.3
        segment n.1: (100, 145)
Domain n.4
        segment n.1: (145, 186)
Domain n.5
        segment n.1: (186, 227)
Domain n.6
        segment n.1: (227, 268)
Domain n.7
        segment n.1: (268, 309)
Domain n.8
        segment n.1: (309, 350)
Domain n.9
        segment n.1: (350, 457)
```