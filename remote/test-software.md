
# haha

> All the install guide assume you are inside the software root directory



for HaMStR

```bash
aria2c http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
aria2c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
aria2c https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
aria2c https://mafft.cbrc.jp/alignment/software/mafft-7.313-gcc_fc6.x86_64.rpm
aria2c -d ~/software https://www.ebi.ac.uk/~birney/wise2/wise2.4.1.tar.gz
```


## HaMStR

`my $prog = 'hmmsearch --cpu 1'`

```bashrc
```

```bash
export WISECONFIGDIR=/usr/share/wise
ln -s /usr/bin/mafft ../bin/mafft-linsi 
cd bin
./configure -n
```


- test

```bash
cd $program/hamstr.v13.2.6/data

../bin/hamstr -h
../bin/

../bin/hamstr -force -sequence_file=testset_cDNA.fa -taxon=test  -hmmset=aves -refspec=Gallus_gallus -hmm=EOG090F0D3R -representative -central
../bin/hamstr -force -sequence_file=testset-prot.fa -taxon=test2 -hmmset=aves -refspec=Gallus_gallus -hmm=EOG090F0D3R -representative -central
../bin/hamstr -force -sequence_file=testset_cDNA.fa -taxon=test  -hmmpath=../core_orthologs -blastpath=../blast_dir -hmmset=modelorganisms_hmmer3 -refspec=DROME -hmm=317.hmm -central  
../bin/hamstr -force -sequence_file=testset_cDNA.fa -taxon=test  -hmmpath=../core_orthologs -blastpath=../blast_dir -hmmset=modelorganisms_hmmer3 -refspec=DROME -hmm=317.hmm -representative -central
../bin/hamstr -force -sequence_file=testset-prot.fa -taxon=test2 -hmmpath=../core_orthologs -blastpath=../blast_dir -hmmset=modelorganisms_hmmer3 -refspec=DROME -hmm=239.hmm -central
../bin/hamstr -force -sequence_file=testset-prot.fa -taxon=test2 -hmmpath=../core_orthologs -blastpath=../blast_dir -hmmset=modelorganisms_hmmer3 -refspec=DROME -hmm=239.hmm -representative -central
../bin/hamstr -force -sequence_file=testset-prot.fa -taxon=test2 -hmmpath=../core_orthologs -blastpath=../blast_dir -hmmset=modelorganisms_hmmer3 -refspec=DROME -hmm=239.hmm -representative -concat -central
```

## [genewise](https://www.ebi.ac.uk/~birney/wise2/)

- test

```
cd $program/wise-2.4.1/test_data

genewise road.pep human.genomic
genewise -hmmer rrm.HMM human.genomic
genewise -hmmer rrm.HMM human.genomic -alg 2193L
estwise road.pep hn_est.fa
estwise -hmmer rrm.HMM hn_est.fa
estwisedb -pfam db.hmm vav.dna
genewisedb -pfam db.hmm vav.dna
```

## augustus

```bash
cd $program/augustus/bin;
./augustus --species=human --UTR=on ../examples/example.fa 
./augustus --species=human --UTR=on ../examples/hsackI10.fa.gz
```

## oases 


```bash
cd $program/oases/velvet/tests
./run-tests.sh

cd $program/oases/data;
oases_pipeline.py -m 21 -M 23 -o singleEnd -d 'test_reads.fa -strand_specific' -p '-min_trans_lgth 100'
cat singleEndMerged/transcripts.fa
```

## Trimmomatic

`software/bin/trimmomatic`
```bash
#!/bin/bash
java -jar $program/Trimmomatic-0.36/trimmomatic-0.36.jar $*
```


## sra-toolkit

[fastq-dump](https://edwards.sdsu.edu/research/fastq-dump/)


