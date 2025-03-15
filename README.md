# find telomeres

`find_telomere.py` is a Python script designed to scan a multiFASTA file for telomere repeats at both the left (start) and right (end) of sequences. The script robustly detects telomeric regions—even when multiple candidate regions exist—by searching for all matches and selecting the best candidate based on the length of the repeat.

## Features

- **Dual-End Scanning:** Checks both the left and right ends of sequences for telomere repeats.
- **Robust Matching:** Uses Python's `re.finditer` to collect all candidate matches and selects the best one based on repeat length and position.
- **Configurable Parameters:** Allows customization of telomere motifs, the minimum number of repeats required, and the search window size.
- **Sorted Output:** Displays results sorted by the number of telomere regions detected (entries with both ends detected appear first).

## Requirements

- Python 3.x

## Installation

Clone or download this repository and navigate to the directory containing the script.


## Usage
Motif does not need to be specified, the default motif is telomeric TTAGGG CCCTAA repeats with min 5 repeats

default:
```shell
./telomere_finder.py input.fasta > telomeres.txt
```

for more customization:
optional arguments
- `--motifs`: Space-separated list of telomere motifs to search for (default: TTAGGG CCCTAA).
- `--min_repeats`: Minimum number of consecutive repeats required (default: 5).
- `--window`: Number of bases from each sequence end to search (default: 200).

```shell
./telomere_finder.py input.fasta --motifs TTAGGG CCCTAA --min_repeats 5 --window 200 > telomeres.txt
```


## Output
By default prints to stdout

Example output:
```
Entry: seq1
  Length: 1500
  Left telomere: YES (TTAGGG) (positions 1-30) sequence: TTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
  Right telomere: YES (CCCTAA) (positions 1471-1500) sequence: CCCTAACCCTAACCCTAACCCTAACCCTAA

Entry: seq2
  Length: 1200
  Left telomere: NO
  Right telomere: YES (TTAGGG) (positions 1171-1200) sequence: TTAGGGTTAGGGTTAGGGTTAGGGTTAGGG

```

## Contact
For questions or suggestions, please open an issue in this repository or contact the maintainer at johannes.nicolaus@gmail.com.
