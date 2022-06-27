# MSAgen
Simple Python module to simulate a set of related DNA sequences

### Dependencies

* [Biopython](https://biopython.org/wiki/Download)
* [Pyvolve](https://github.com/sjspielman/pyvolve)

### Usage

Example:
```python
import MSAgen

sequences = MSAgen.generate_sequences(num_sequences = 10, # the number of sequences to generate 
                                      seqlen = 10000, # length of each sequence (in bp)
                                      genelen = 2000, # length of the gene in each sequence (in bp, can be 0)
                                      coding_dist = 0.2, # branch length of the underlying tree for simulated gene evolution
                                      noncoding_dist = 0.4) # branch length for flanking regions
```

### Limitations

Currently, only no or one single-exon gene can be generated, possibly flanked by less conserved regions (flanked regions are created when `seqlen` > `genelen`, the genes are then centered in the sequences)
