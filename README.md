A pipeline for BCR repertoire libraries from  - UMI Barcoded Illumina MiSeq 325+275 paired-end 5’RACE BCR mRNA.


Library preperation and sequencing method:

The sequences were amplified specific primers 1.isotypes human Primers 2. SmartNNNa Primers. 
The generated libraries were then sequenced with Illumins MiSeq 325+275.

Input files:

* Pair-end reads Sample_R1.fastq and Sample_R2.fastq - The input files are from the `Yaari1` module, which process illumina raw reads into samples pair-end fastq, based on the M1S and Z pipeline. 
* primers sequences - SmartNNNa.fasta , isotypes.human.fasta.
* Assemble pairs reference

Output file:

1. Sample.fasta - Single processed fast file
2. log tab file for each steps
3. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

1. MaskPrimer align
2. PairAwk
3. FilterSeq quality
4. FilterSeq length
4. Cluster UMIs
5. AlignSets muscle
6. BuildConsensus
7. PairSeq
8. AssemblePairs sequential


Primers used:

* [SmartNNNa Primers (R1)](add link to git)]

* [isotypes human Primers (R2)](add link to git)



Refrence used:

*[humanIGHV w gaps](...)