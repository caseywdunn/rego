# rego: REgular expression extraction of sequence reads with oliGO

A tool to extract sequence reads based on oligo sequence (like primers) using regular expressions.

## Example

```bash
python rego.py --output Agalma_elegans.18s TGGTTGATCCTGCCAGT data/Agalma_elegans.18s.64.fastq
```


## Test data

`Agalma_elegans.18s.64.fastq` - 64 sequence reads from the siphonophore *Agalma elegans* that contain the 18S sequence `TGGTTGATCCTGCCAGT`

## Development

To run tests:

    python rego.py GG null.fastq --test

The first two arguments are ignored when running a test and can be any strings.