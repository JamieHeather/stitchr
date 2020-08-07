#stiTChR 0.4.1

### Stitch together TCR coding nucleotide sequences from V/J/CDR3 information

##### Jamie Heather | Cobbold lab @ MGH | 2018

### Summary

Sometimes you need a TCR nucleotide sequence, but all you have is limited information. This script aims to generate *a*\* coding nucleotide sequence for a given rearrangement (e.g. for use when generating TCR expression vectors) in those situations.

The script takes the known V/J/CDR3 information, and uses that to pull out the relevant germline TCR nucleotide sequences and stitch them together, using common human codons to fill in the non-templated regions.

### Installation and dependencies

As of version 0.4 ```stiTChR``` is designed to be run on Python 3.8, and has primarily been tested on Python 3.7.7 and 3.8.3. 

Simply clone the repo to a desired location, navigate to the Scripts directory, then you can run the script via the command line as detailed below.

The only non-default package used is Biopython, which can be easily installed via `pip`:

```bash
pip3 install biopython
```

## Example usage 

```stiTChR``` uses relative paths. Please ensure you are in the Scripts directory to run the script. 

```bash
python3 stitchr.py -v [IMGT V gene] -j [IMGT J gene] -cdr3 [CDR3aa]

python3 stitchr.py -v TRBV7-3*01 -j TRBJ1-1*01 -cdr3 CASSYLQAQYTEAFF

python3 stitchr.py -v TRAV1-2 -j TRAJ33 -cdr3 CAVLDSNYQLIW
```

### Usage notes

*This script will almost certainly not produce the original recombined sequences that encoded the original TCRs; it merely recreates an equivalent full length sequence that will encode the same amino acid sequences. It aims to produce a sequence as close to germline as possible - so all CDR3 residues that *can* be germline encoded are such, in this ouput. Non-templated residues in the CDR3 (or those templated by the D gene, which is treated as non-templated here) are chosen from known codon usage frequencies.

This script currently only works on human alpha/beta TCR chains, but should be readily adapted to any other species/locus.

Care must be taken to ensure that the correct TCR informaton is input, e.g. ensure you're using proper IMGT gene nomenclature (older deprecated gene names will not work), and have the correct and full CDR3 sequence from the conserved cysteine to the conserved phenylalanine (or rarely, tryptophan) residue. The script produces a TCR from the information given, trying to provide warnings or errors if it detects an improbable or implausible combination, yet it's possible that the script might produce output that *looks* OK yet which does not reproduce a coding sequence for the intended TCR. 

If you request an allele for which there isn't complete sequence data, the script will attempt to default to the prototypical allele (*01) of that gene. If it cannot find sequence for that then it will throw an error. Similarly it will attempt to use the correct leader seqeunces (L-PART1+L-PART2) for the specified allele, but if it can't find one it'll default back to the prototype's. Note that IMGT-provided gene sequences which are 'partial' at either end of their sequence are discounted entirely, as full length sequences are needed.

By default the script will use the TRBC gene located in the same cluster as the J gene (i.e. TRBJ1-1 through TRBJ1-6 will get TRBC1, while TRBJ2-1 through TRBJ2-7 will get TRBC2). This can be overriden (see optional arguments).

All required files are included in the repo. If you want to change or update the IMGT data files, you'll need to re-run `split-imgt-data.py`.

```StiTChR``` can be run in a higher-throughput mode, using a tab-separated input file - see the instructions for [**```thimble```**](#Thimble) below.


### Optional arguments

* `-h` - see a help menu, containing all the command line options
* `-c` - specify a particular constant region gene (in the case of TRBC) or allele
* `-s` - specify a species: 'human' or 'mouse' are the only valid options currently, with human as default 
* `-aa` - provide an incomplete amino acid sequence (spanning at least the CDR3, with some padding on either side), to assess the accuracy of the stitched TCR sequence. Must be a single string, unbroken by spaces or linebreaks
* `-cu` - use an alternative codon usage file, from which to generate the sequences for the non-templated residues (see below)
* `-l` - use a different leader region to that present with the given V  
* `-n` - provide a name for the TCR chain, which will be included in the FASTA file header

#### Codon usage files

Non-templated based are assigned by taking the most common nucleotide triplet for a given amino acid, in a provided codon usage file.

The default codon usage files are taken straight from the default Kazusa [human](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606) (Homo sapiens [gbpri]: 93487) and [mouse](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10090) (Mus musculus [gbrod]: 53036) entries. Alternative files can be provided, but must be in the same format (e.g. those provided by [HIVE](https://hive.biochemistry.gwu.edu/dna.cgi?cmd=refseq_processor&id=569942)). U/T can be used interchangeably, as all U will be replaced with T anyway.

#### Providing a partial amino acid sequence

If you provide a partial amino acid sequence ```stiTChR``` will perform a rudimentary pairwise alignment, just to give a quick visual assessment of the quality of the sequence generation.

#### A fancier example

Let's take the example of the well described A2-NLV restricted [C25 TCR from the 5D2N PDB structure](https://www.rcsb.org/structure/5d2n). We can take the amino acid sequence straight from the PDB FASTA file:

```
>5D2N:E|PDBID|CHAIN|SEQUENCE
MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTI
QRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVEL
SWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSA
EAWGRAD
```

We can then pull out the V, J, and CDR3 information. There's lots of ways to do this, but the easiest manual way is to find the CDR3 and then search the immediately neighbouring sequences against V/J amino acid sequences (obtainable via IMGT/GENE-DB). This gives:
  
TRBV7-6 / TRBJ1-4 / CASSLAPGTTNEKLFF

Then we can run the code like this:

```bash
python3 stitchr.py -v TRBV7-6 -j TRBJ1-4 -cdr3 CASSLAPGTTNEKLFF -n C25 -aa MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

# Produces the following output

>nt-C25-TRBV7-6*01-TRBJ1-4*01-TRBC1*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
ATGGGCACCAGTCTCCTATGCTGGGTGGTCCTGGGTTTCCTAGGGACAGATCACACAGGTGCTGGAGTCTCCCAGTCTCCCAGGTACAAAGTCACAAAGAGGGGACAGGATGTAGCTCTCAGGTGTGATCCAATTTCGGGTCATGTATCCCTTTATTGGTACCGACAGGCCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCAATTATGAAGCCCAACAAGACAAATCAGGGCTGCCCAATGATCGGTTCTCTGCAGAGAGGCCTGAGGGATCCATCTCCACTCTGACGATCCAGCGCACAGAGCAGCGGGACTCGGCCATGTATCGCTGTGCCAGCAGCCTGGCCCCCGGCACCACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC

>aa-C25-TRBV7-6*01-TRBJ1-4*01-TRBC1*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF

MG------------------AGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ
||                  ||||||||||||||||||||||||||||||||||||||||
MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ

GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT

NEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFY-PDHVELS
|||||||||||||||||||||||||||||||||||||||||||||||||||  |||||||
NEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGF-FPDHVELS

WWVNGKEVHSGVC-TDPQPLKEQPALNDSRYA-LSSRLRVSATFWQNPRNHFRCQVQFYG
||||||||||||  |||||||||||||||||  |||||||||||||||||||||||||||
WWVNGKEVHSGV-STDPQPLKEQPALNDSRY-CLSSRLRVSATFWQNPRNHFRCQVQFYG

LSENDEWTQDRAKPVTQIVSAEAWGRAD--------------------------------
||||||||||||||||||||||||||||                                
LSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVL

---------------
               
VSALVLMAMVKRKDF
```

We can see that there's a few mismatches in the latter half of the stitched sequence, so perhaps this crystal actually used the other TRBC gene. We can swap that in:

```bash
python3 stitchr.py -v TRBV7-6 -j TRBJ1-4 -cdr3 CASSLAPGTTNEKLFF -n C25 -c TRBC2 -aa MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

# Produces:

>nt-C25-TRBV7-6*01-TRBJ1-4*01-TRBC2*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
ATGGGCACCAGTCTCCTATGCTGGGTGGTCCTGGGTTTCCTAGGGACAGATCACACAGGTGCTGGAGTCTCCCAGTCTCCCAGGTACAAAGTCACAAAGAGGGGACAGGATGTAGCTCTCAGGTGTGATCCAATTTCGGGTCATGTATCCCTTTATTGGTACCGACAGGCCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCAATTATGAAGCCCAACAAGACAAATCAGGGCTGCCCAATGATCGGTTCTCTGCAGAGAGGCCTGAGGGATCCATCTCCACTCTGACGATCCAGCGCACAGAGCAGCGGGACTCGGCCATGTATCGCTGTGCCAGCAGCCTGGCCCCCGGCACCACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTCACCTCCGAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTCAAGAGAAAGGATTCCAGAGGC

>aa-C25-TRBV7-6*01-TRBJ1-4*01-TRBC2*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG

MG------------------AGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ
||                  ||||||||||||||||||||||||||||||||||||||||
MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ

GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT

NEKLFFGSGTQLSVLEDLNK-VFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELS
|||||||||||||||||| | |||||||||||||||||||||||||||||||||||||||
NEKLFFGSGTQLSVLEDL-KNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELS

WWVNGKEVHSGVC-TDPQPLKEQPALNDSRYA-LSSRLRVSATFWQNPRNHFRCQVQFYG
||||||||||||  |||||||||||||||||  |||||||||||||||||||||||||||
WWVNGKEVHSGV-STDPQPLKEQPALNDSRY-CLSSRLRVSATFWQNPRNHFRCQVQFYG

LSENDEWTQDRAKPVTQIVSAEAWGRAD--------------------------------
||||||||||||||||||||||||||||                                
LSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVL

-----------------
                 
VSALVLMAMVKRKDSRG
```

This produces even more mismatches - so the constant region used in the crystal has presumably been altered for expression/crystallization purposes.

#### A mouse example

```StiTChR``` also now supports murine TCRs, however it should be noted that due to poorer quality annotations in IMGT the sequences produced should be treated with even more caution than used for human data.

Here's an example of how to run ```stiTChR``` on everyone's favourite mouse TCR, OT-I (with the actual sequence inferred from [this plasmid on AddGene](https://www.addgene.org/52111/):

```bash
python3 stitchr.py -s mouse -v TRBV12-1 -j TRBJ2-7 -cdr3 CASSRANYEQYF
python3 stitchr.py -s mouse -v TRAV14-1 -j TRAJ33 -cdr3 CAASDNYQLIW 
```

#### A note on CDR3 C-terminal residues

```StiTChR``` assumes that the J gene will not undergo deletion past the C-terminal residue of the CDR3 (which occurs approximately in the middle of the J). Thus the code looks for the appropriate residue at the end of the CDR3, which in the majority of cases will be a phenylalanine (F). However in some cases it might be something else, like a W (not uncommon in human TRAJ/mice genes) or even something more exotic like a C, L or J (which occur in certain mouse J genes). Note that most of these non-F/W residues are found in J genes with a predicted ['ORF' IMGT status](http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html), and thus might not contribute to functioning TCRs, but stiTChR will still let you generate a plausible sequence using them.

# Thimble  0.1.0

### Run stiTChR high-throughput on multiple and paired TCRs

Instead of running ```stiTChR``` on rearrangements one by one, you can fill out the necessary details into a tab separated file (.tsv) and submit it to ```thimble```.

The format of the input data is specified in the file 'bulk_input_template.tsv', located in the root directory, with some examples shown in 'bulk_input_example.tsv'. All of the recombination-specific fields that can ordinarily be specified at the command line in ```stiTChR``` can be applied per row using ```thimble```, with the exception of species which must be kept the same for all TCRs in a given run.

Note that the input to ```thimble``` can also be used to generate rearrangements for both the alpha and beta chain of a given clonotype on one row, with additional options to link those sequences together. A number of x2A potential linkers are provided (in the Data/linkers.tsv file). If custom linkers are desired, you can either edit that linkers file or just enter the nucleotide sequence of the desired linker into the Linker column of the input tsv.

Any warnings and errors generated on a per-TCR basis are recorded in the final output file; it is recommended that users check this information, to ensure they understand the limitations of a specific sequence.  

## Example usage 

```bash
python3 thimble.py -in [input tsv] -o [output tsv] 

python3 thimble.py -in ../bulk_input_example.tsv -o testing 
```

### Optional arguments

* `-h` - see a help menu, containing all the command line options
* `-s` - specify a species: 'human' or 'mouse' are the only valid options currently, with human as default 
* `-cu` - use an alternative codon usage file, from which to generate the sequences for the non-templated residues (see below)

