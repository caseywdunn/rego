import argparse
import re
import unittest
import numpy as np
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
    
    def get_perfect_overlap(self, other):
        # Given another ReadHit, return the length of perfect overlap between the two
        # The match must be perfect over the entire length of the overlap 
        # Consider every possible overlap, and return 0 if no perfect overlap is found

        overlap = 0
        # Check if one read is a subset of the other
        if self.read in other.read or other.read in self.read:
            overlap = min(len(self.read), len(other.read))
        else:
            for i in range(0, len(self.read)):
                if self.read[i:] == other.read[:len(self.read) - i]:
                    overlap = max(overlap, len(self.read) - i)

            for i in range(0, len(other.read)):
                if other.read[i:] == self.read[:len(other.read) - i]:
                    overlap = max(overlap, len(other.read) - i)
        return overlap

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
    
    def test_get_perfect_overlap(self):
        read1 = ReadHit("AGTAACAAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCG", 0, 8)
        read2 = ReadHit(      "AAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCGAGTGCC", 0, 7)
        self.assertEqual(read1.get_perfect_overlap(read2), 42)
        self.assertEqual(read2.get_perfect_overlap(read1), 42)

        read1 = ReadHit("AGTAACAAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCG", 0, 8)
        read2 = ReadHit(      "AAACGCGTGTAC", 0, 7)
        self.assertEqual(read1.get_perfect_overlap(read2), 12)
        self.assertEqual(read2.get_perfect_overlap(read1), 12)

        read1 = ReadHit("AGTAACAAACGCGTGTACT", 0, 8)
        read2 = ReadHit("GGGATTCGATTCCCTTAACAGGTCAGTCG", 0, 7)
        self.assertEqual(read1.get_perfect_overlap(read2), 0)
        self.assertEqual(read2.get_perfect_overlap(read1), 0)

    def test_find_components(self):
        hits = [ReadHit("AGTAACAAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCG", 0, 8),
                ReadHit(      "AAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCGAGTGCC", 0, 7),
                ReadHit("GGGGGGGGGCCCCCCCTTTTTTTTAAAAAA", 0, 7)]
        self.assertEqual(find_components(hits, 20), [[0, 1], [2]])
    
    def test_get_consensus(self):
        lines = ["AGTAACAAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCG----------",
                 "------AAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCGAGTGCC----",
                 "------------GTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCGAGTGCCAATG"]
        self.assertEqual(get_consensus(lines), "AGTAACAAACGCGTGTACTGGGATTCGATTCCCTTAACAGGTCAGTCGAGTGCCAATG")




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

def dfs(node, adj_matrix, visited, component):
    visited[node] = True
    component.append(node)
    for neighbor, is_connected in enumerate(adj_matrix[node]):
        if is_connected and not visited[neighbor]:
            dfs(neighbor, adj_matrix, visited, component)

def find_connected_components(adj_matrix):
    n = len(adj_matrix)  # Number of nodes
    visited = [False] * n
    components = []
    
    for node in range(n):
        if not visited[node]:
            component = []
            dfs(node, adj_matrix, visited, component)
            components.append(component)
    
    return components

def find_components(hits, overlap):
    # create adjacency matrix
    n = len(hits)
    adj = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            adj[i, j] = hits[i].get_perfect_overlap(hits[j])

    # regularize adjacency matrix
    adj = np.where(adj >= overlap, 1, 0)

    # find connected components
    components = find_connected_components(adj)
    return components

def get_consensus(lines):
    # Create a numpy array where each row is a line
    # and each column is a nucleotide
    # The lines are already padded

    # Assert that all lines are the same length and nonzero length
    line_length = len(lines[0])
    assert line_length > 0
    for line in lines:
        assert len(line) == line_length

    n = len(lines)
    m = len(lines[0])
    matrix = np.zeros((n, m), dtype=np.uint8)
    for i in range(n):
        for j in range(m):
            matrix[i, j] = ord(lines[i][j])

    # Create a consensus sequence
    consensus = ""
    for j in range(m):
        counts = dict()
        for i in range(n):
            if matrix[i, j] in counts:
                counts[matrix[i, j]] += 1
            else:
                counts[matrix[i, j]] = 1
        
         # ignore - unless it is the only character
        max_count = 0
        max_nuc = 0
        counts.pop(ord('-'), None)
        if len(counts) == 0:
            consensus += '-'
        else:
            for nuc, count in counts.items():
                if count > max_count:
                    max_count = count
                    max_nuc = nuc
            consensus += chr(max_nuc)

    return consensus

def main(oligo, input_file, output_file, before, after, max_hits):
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

    # Find components and create a list of lists of hits that are connected
    components = find_components(hits, 20)

    hits_grouped = []
    for component in components:
        hits_grouped.append([hits[i] for i in component])

    consensuses = dict()
    
    for c, group in enumerate(hits_grouped):
        print(f"Component {c}")
        for hit in group:
            print(hit.get_seq())
        print()


        # Loop through the hits and get the maximum number of characters before and after oligo_start
        max_before = 0
        max_after = 0
        for hit in group:
            max_before = max(max_before, hit.oligo_start)
            max_after = max(max_after, len(hit.read) - hit.oligo_start)
        
        # create hit_lines of output text
        hit_lines = []
        longest_line = 0
        for hit in group:
            new_line = f"{'-' * (max_before - hit.oligo_start)}{hit.get_seq()}"
            longest_line = max(longest_line, len(new_line))
            hit_lines.append(new_line)
        
        # Loop through the hit_lines and
        # pad the trailing ends with - to make them all the same length
        for i in range(len(hit_lines)):
            hit_lines[i] = hit_lines[i].ljust(longest_line, '-')

        # Create a string that has all the lines in phylip format
        phylip_text = f"{len(hit_lines)} {longest_line}\n"
        for i in range(len(hit_lines)):
            label = f"read_{i}".ljust(12)
            phylip_text += f"{label}  {hit_lines[i]}\n"
        print(f"Phylip format for component {c}")
        print(phylip_text)

        # Write the output to a file
        if output_file:
            with open(output_file + f".component_{c}.reads.phy", 'w') as f:
                f.write(phylip_text)
        
        consensus = get_consensus(hit_lines)
        consensuses[c] = consensus
    
    # Write the consensus to a fasta file
    if output_file:
        with open(output_file + ".consensus.fasta", 'w') as f:
            for c, consensus in consensuses.items():
                f.write(f">component_{c}\n")
                f.write(consensus + "\n")
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Use regular expressions with an oligo to extract a specified region from a sequence. The oligo can contain IUPAC degenerate nucleotide codes.")
    parser.add_argument("oligo", type=str, help="the oligo sequence")
    parser.add_argument("input", type=str, help="input fastq file")
    parser.add_argument("--output", type=str, help="Name to prepend to output, can include a path", default=None)
    parser.add_argument("--before", type=int, help="number of nucleotides to preserve before the oligo", default=None)
    parser.add_argument("--after", type=int, help="number of nucleotides to preserve after the oligo", default=None)
    parser.add_argument("--max_hits", type=int, help="maximum number of read hits to retrieve", default=20)
    parser.add_argument("--test", action="store_true", help="run tests")
    
    args = parser.parse_args()
    
    if args.test:
        unittest.main(argv=[''], exit=False)
        exit()

    main(args.oligo, args.input, args.output, args.before, args.after, args.max_hits)