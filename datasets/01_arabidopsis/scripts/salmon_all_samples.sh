# get salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz
tar xf SalmonBeta-0.6.0_DebianSqueeze.tar.gz
mv SalmonBeta-0.6.1_DebianSqueeze/* .

# get reference
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-30/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.30.cdna.all.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.30.cdna.all.fa.gz

# index reference
./bin/salmon index -t Arabidopsis_thaliana.TAIR10.30.cdna.all.fa -i transcripts_index --type quasi -k 31

# run quant
SAMPLES=(SRR3136731 SRR3136732 SRR3136733 SRR3136734 SRR3136735 SRR3136736 SRR3136737 SRR3136738 SRR3136739)
for id in ${SAMPLES[@]}; do
  ./bin/salmon quant \
    -p 20 \
    -i transcripts_index \
    -l U \
    -r <(gunzip -c ${id}.fastq.gz) \
    -o quant_${id}
; done

# bundle up the results
tar zcf salmon_output.tar.gz quant_*
