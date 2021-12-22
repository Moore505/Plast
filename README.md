# Plast (Primitive Local Alignment Search Tool)
## December 6th 2021
## IFT3295 - Bioinformatics
## Université de Montréal

### Description of plast
Plast is a heuristic algorithm freely based on the original, "gap-free" version of blast for nucleotide sequences. More details on Blast [here](https://en.wikipedia.org/wiki/BLAST_(biotechnology)).

### Packages required
- **Numpy**

### Tested versions of Python
- **3.9.7**

### Exemple of using plast with default arguments
To use Plast with the default optionnal arguments, the program must be called with the database file name and the query to search, here's an example:

``` python
python plast.py -db "tRNAs.fasta" -i "CGTAGTCGGCTAACGCATACGCTTGATAAGCGTAAGAGCCC"
```

Here are the details for the arguments of Plast :

| Argument | Requirement | Default value | Description
| --- | ----------- | ----------- | ----------- |
| -db | Mandatory | - | file name of the database of known sequences 
| -i | Mandatory | - | query
| -E | Optionnal | 4 | treshold
| -ss | Optionnal | 0.001 | Signification threshold
| -seed | Optionnal | '11111111111' | Seed for the search (non-functionnal right now))