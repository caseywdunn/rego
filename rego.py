import argparse
import re
import unittest
# A contraction of "regular expression with oligo"
# written by: Casey Dunn, casey.dunn@yale.edu, https://dunnlab.org

degenerate_nucleotides = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]'
}

def primer_to_regex(primer):
    regex = ''.join(degenerate_nucleotides[nuc] for nuc in primer)
    return regex

def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'H': 'D', 'D': 'H', 'N': 'N'}
    return ''.join(complement[nuc] for nuc in reversed(sequence))

class ReadHit():
    def __init__(self, read, oligo_start, oligo_end):
        self.read = read
        self.oligo_start = oligo_start
        self.oligo_end = oligo_end
    
    def get_seq(self, lowercase=True):
        if lowercase:
            read = self.read.lower()
            # Make the oligo uppercase
            read = read[:self.oligo_start] + read[self.oligo_start:self.oligo_end].upper() + read[self.oligo_end:]
            return read
        else:
            return self.read.upper()
        return self.read[self.start:self.end]


def merge (seq1, seq2, overlap=15):
    # Given two sequences, return the merged sequence if seq2 extends seq1
    # and overlaps the end of seq1 by at least the specified number of nucleotides
    # Otherwise, return an empty string
    for i in range(0, len(seq1) - overlap + 1):
        if seq2.startswith(seq1[i:]):
            return seq1[:i] + seq2
    return ''
    
def overlap_assembler(sequences, overlap=15):
    # Given a list of sequences, return a list of merged sequence
    # Stop when no more sequences can be merged
    

    while True:
        merged = False
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i != j:
                    merged_seq = merge(sequences[i], sequences[j], overlap)
                    if merged_seq:
                        sequences[i] = merged_seq
                        sequences.pop(j)
                        merged = True
                        break
            if merged:
                break
        if not merged:
            break
    return sequences

class TestRegoMethods(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertEqual(reverse_complement("ATCGR"), "YCGAT")
        #self.assertEqual(1,2)

    def test_merge(self):
        seq1 = "ATCGTACG"
        seq2 = "TACGCGT"
        self.assertEqual(merge(seq1, seq2, 3), "ATCGTACGCGT")
        self.assertEqual(merge(seq2, seq1, 3), "")
        self.assertEqual(merge("CGTAGGTA", "AAAAAAAAA", 3), "")
    
    def test_overlap_assembler(self):
        sequences = ["ATCGTACG", "TACGCGT", "CGTACGCGT"]
        self.assertEqual(overlap_assembler(sequences, 3), ["ATCGTACGCGT"])
        self.assertEqual(overlap_assembler(["CGTAGGTA", "AAAAAAAAA", "GGGGGGGGG"], 5), ["CGTAGGTA", "AAAAAAAAA", "GGGGGGGGG"])


def trim_read(oligo, read, before, after):
    forward_regex = primer_to_regex(oligo)
    
    forward_match = re.search(forward_regex, read)

    read = read.lower()
    
    if forward_match:
        start = forward_match.start()
        end = forward_match.end()
        if before:
            start = max(0, start - before)
        if after:
            end = min(len(read), end + after)

        # Make the oligo uppercase
        read = read[:start] + read[start:end].upper() + read[end:]

        return read[start:end]
    else:
        return read

def main(oligo, input_file, before, after, max_hits):
    forward_regex = primer_to_regex(oligo)
    reverse_regex = primer_to_regex(reverse_complement(oligo))
    
    hits = []
    i = 0
    with open(input_file, 'r') as f:
        for line in f:
            i += 1
            if i % 4 != 2:
                continue
            
            if len(hits) >= max_hits:
                break
            
            if re.search(reverse_regex, line):
                line = line.strip()
                line = reverse_complement(line)

            forward_re = re.search(forward_regex, line)
            if forward_re:
                read = line.strip()
                oligo_start = forward_re.start()
                oligo_end = forward_re.end()
                hit = ReadHit(read, oligo_start, oligo_end)
                hits.append(hit)
                
    print(f"Hits found: {len(hits)}")

    # Sort the hits by hit.oligo_start
    hits.sort(key=lambda x: x.oligo_start)

    # Loop through the hits and get the maximum number of characters before and after oligo_start
    max_before = 0
    max_after = 0
    for hit in hits:
        max_before = max(max_before, hit.oligo_start)
        max_after = max(max_after, len(hit.read) - hit.oligo_start)
    
    # Print the hits, but pad with spaces so that the oligo is in the same position
    for hit in hits:
        print(f"{' ' * (max_before - hit.oligo_start)}{hit.get_seq()}")
    print("\n")
        
    assembled = overlap_assembler(hits)
    assembled.sort(key=len, reverse=True)
    print("Assembled sequences:")
    
    n = 0
    for sequence in assembled:
        print(f"> Sequence {n}:")
        print(sequence)
        n += 1
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Use regular expressions with an oligo to extract a specified region from a sequence. The oligo can contain IUPAC degenerate nucleotide codes.")
    parser.add_argument("oligo", type=str, help="the oligo sequence")
    parser.add_argument("input", type=str, help="input fastq file")
    parser.add_argument("--before", type=int, help="number of nucleotides to preserve before the oligo", default=None)
    parser.add_argument("--after", type=int, help="number of nucleotides to preserve after the oligo", default=None)
    parser.add_argument("--max_hits", type=int, help="maximum number of read hits to retrieve", default=20)
    parser.add_argument("--test", action="store_true", help="run tests")
    
    args = parser.parse_args()
    
    if args.test:
        unittest.main(argv=[''], exit=False)
        exit()

    main(args.oligo, args.input, args.before, args.after, args.max_hits)