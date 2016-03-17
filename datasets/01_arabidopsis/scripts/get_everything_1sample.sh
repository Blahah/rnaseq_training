# get Salmon and unpack it
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz
tar xf SalmonBeta-0.6.0_DebianSqueeze.tar.gz
mv SalmonBeta-0.6.1_DebianSqueeze/* .
rm -rf SalmonBeta*

# get just the first read set
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/002/SRR3136732/SRR3136732.fastq.gz

# get reference
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-30/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.30.cdna.all.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.30.cdna.all.fa.gz
