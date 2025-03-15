#!/usr/bin/env python3
import argparse
import re

def find_telomere_end(seq, motifs, min_repeats, window, end='left', max_offset=10):
    """
    Look for telomere motifs in a specific sequence end.
    This updated version collects all candidate matches and selects the best one based on the length of the repeat.
    
    For the left end:
      - A candidate is valid if it starts within 'max_offset' bases from the start.
      
    For the right end:
      - A candidate is valid if its end is within 'max_offset' bases of the end of the search window.
      
    Returns:
      found (bool), match sequence (str), motif (str), start (int), end (int)
    """
    candidates = []
    
    if end == 'left':
        region = seq[:window]
        for motif in motifs:
            pattern = f"(?:{motif}){{{min_repeats},}}"
            for match in re.finditer(pattern, region):
                if match.start() <= max_offset:
                    candidates.append((match, motif))
        if candidates:
            # Choose the candidate with the longest match; if tied, the one that starts earliest.
            best = max(candidates, key=lambda x: (len(x[0].group()), -x[0].start()))
            match, motif = best
            return True, match.group(), motif, match.start(), match.end()
    
    elif end == 'right':
        region = seq[-window:]
        for motif in motifs:
            pattern = f"(?:{motif}){{{min_repeats},}}"
            for match in re.finditer(pattern, region):
                # Valid if the match ends close to the end of the region.
                if (len(region) - match.end()) <= max_offset:
                    candidates.append((match, motif))
        if candidates:
            best = max(candidates, key=lambda x: (len(x[0].group()), x[0].end()))
            match, motif = best
            offset = len(seq) - window
            return True, match.group(), motif, match.start() + offset, match.end() + offset

    return False, None, None, None, None

def process_entry(seq_id, seq, motifs, min_repeats, window):
    left_found, left_match, left_motif, left_start, left_end = find_telomere_end(seq, motifs, min_repeats, window, end='left')
    right_found, right_match, right_motif, right_start, right_end = find_telomere_end(seq, motifs, min_repeats, window, end='right')

    # Score: count how many ends have telomere repeats (2, 1, or 0)
    score = int(left_found) + int(right_found)
    
    return {
        "seq_id": seq_id,
        "length": len(seq),
        "left_found": left_found,
        "left_match": left_match,
        "left_motif": left_motif,
        "left_start": left_start,
        "left_end": left_end,
        "right_found": right_found,
        "right_match": right_match,
        "right_motif": right_motif,
        "right_start": right_start,
        "right_end": right_end,
        "score": score
    }

def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id is not None:
            yield seq_id, "".join(seq_lines)

def main():
    parser = argparse.ArgumentParser(
        description="Scan a FASTA file for telomere repeats at sequence ends and sort entries by telomere presence."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("--motifs", nargs="+", default=["TTAGGG", "CCCTAA"],
                        help="Telomere motifs (default: TTAGGG CCCTAA)")
    parser.add_argument("--min_repeats", type=int, default=5,
                        help="Minimum consecutive repeats required (default: 5)")
    parser.add_argument("--window", type=int, default=200,
                        help="Window size (in bp) at each end to search (default: 200)")
    args = parser.parse_args()

    results = []
    for seq_id, seq in parse_fasta(args.fasta):
        result = process_entry(seq_id, seq, args.motifs, args.min_repeats, args.window)
        results.append(result)

    # Sort entries by the "score" (number of ends with telomere repeats) in descending order.
    results = sorted(results, key=lambda r: r["score"], reverse=True)

    for res in results:
        print(f"Entry: {res['seq_id']}")
        print(f"  Length: {res['length']}")
        if res["left_found"]:
            print(f"  Left telomere: YES ({res['left_motif']}) (positions {res['left_start']+1}-{res['left_end']}) sequence: {res['left_match']}")
        else:
            print("  Left telomere: NO")
        if res["right_found"]:
            print(f"  Right telomere: YES ({res['right_motif']}) (positions {res['right_start']+1}-{res['right_end']}) sequence: {res['right_match']}")
        else:
            print("  Right telomere: NO")
        print("")

if __name__ == "__main__":
    main()
            
