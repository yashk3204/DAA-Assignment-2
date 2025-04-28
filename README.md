# DAA Assignment - 2
## Group - 32

    1. Yash Kantamneni (2022A7PS0120H)
    2. Anshul Sharma (2022A7PS1290H)
    3. Raghav Singh Sengar (2022A7PS1797H)
    4. Sahil Yelventage (2022A7PS1351H)
    5. Aryan Deshmukh (2022A7PS1283H)

## Execution Instructions:

`1.  EXACT (ALGORITHM 1)`
- Keep the datasets in the same directory as the code "exact.cpp", and change the path of the infile in the main function (line 636) to this .txt file.
- Also change the path to the outfile to any .txt file in line 642 of the exact.cpp file.
- Run it using gcc/g++ and the output will be printed in a new .txt file given by the user in line 642.

`2. CORE EXACT (ALGORITHM 4)`
- Keep the datasets in the same directory as the code "coreExact.cpp", and change the path of the infile in the main function (line 300) to this .txt file.
- Run it using gcc/g++ and the output will be printed on the terminal.

## Dataset Preparation

- The datasets were preprocessed by sorting all unique node IDs and then re-indexing them using a map.
- We also removed edges of the type (a,a) i.e. self-loops.
- We also removed dupicate undirected edges i.e. if (a,b) and (b,a) are both present, we removed one of them.
- We stored all the edges in a set.
- To our preprocessed .txt file, we first printed the number of nodes and number of edges in the first line, followed by all the edges consisting of re-indexed vertices (1-indexed).
- Refer to the "preprocessor.py" file in the submitted folder or in the GitHub link for the code.

## Links

- [Report](https://drive.google.com/file/d/1F-Liryw913W45d5psYy_tRn-V_4UHM0L/view?usp=sharing)
- [GitHub Repository](https://github.com/yashk3204/DAA-Assignment-2)
- [Drive link to datasets - raw and processed](https://drive.google.com/drive/folders/1viBjIm0UYmm124JUhJw6uPYTB3F0TBuf?usp=sharing)

## Results

### Algorithm 1 - Exact
| Dataset | H-Value | No. of (h-1) Cliques | No. of Triangles | No. of Edges in Densest Subgraph | No. of Vertices in Densest Subgraph | Density | Runtime |
| :------ | :-------| :--------------------| :--------------- | :------------------------------- | :---------------------------------- | :------ | :------ |
| AS-733 | 3 | 3132 | 2503 | 226 | 28 | 8.0714 | 61.3549s |
| CA-HepTh | 3 | 25973 | 28339 | 496 | 32 | 15.5 | 619.24s |
| Yeast | 3 | 1948 | 206 | 19 | 7 | 2.7142 | 5.26113s |

### Algorithm 4 - CoreExact
| Dataset | H-Value | No. of (h-1) Cliques | No. of Triangles | No. of Edges in Densest Subgraph | No. of Vertices in Densest Subgraph | Density | Runtime |
| :-------| :-------| :--------------------| :--------------- | :------------------------------- | :---------------------------------- | :------ | :------ |
| AS-733 | 3 | 3132 | 2503 | 556 | 168 | 3.3095 | 0.4194s |
| CA-HepTh | 3 | 25973 | 28339 | 22842 | 6882 | 3.3191 | 2.47503s |
| Yeast | 3 | 1948 | 206 | 406 | 236 | 1.7203 | 0.07966s |
| AS-Caida | 3 | 52861 | 34265 | 26549 | 8393 | 3.1632 | 32.0799s |

## Individual Contributions

| Name | Contribution |
| :-------- | :------- |
| Yash Kantamneni | Code for algorithm-1 (EXACT) and its report/observations, GitHub repo |
| Anshul Sharma | Data preparation/preprocessing and project webpage |
| Raghav Singh Sengar | Code for algorithm-4 (CORE-EXACT) and its report/observations, project webpage |
| Sahil Yelvantge | Code for algorithm-1 (EXACT) and its report/observations, and README |
| Aryan Deshmukh | Code for algorithm-4 (CORE-EXACT) and its report/observations, and README |