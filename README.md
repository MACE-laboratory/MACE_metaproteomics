# Metaproteomics pipeline for the MACE laboratory
Massimo Bourquin, 2024

# 1. Raw to mzML conversion
To use SAGE for DIA LFQ: 
```
thermorawfileparser -i SAMPLE.raw  -N -f 1 -g
```
This results into the mzML.gz files that can be used to run sage.

To use DIA-NN for DIA LFQ: 
```
mscli SAMPLE.raw -o output_folder --mzML --64 --zlib --filter "peakPicking vendor msLevel=1-"
```
This results into the mzML.gz files that can be used to run DIA-NN.

# 2.a Running Sage to perform LFQ
example: 
```
sage config.json --batch-size 12
```

# 2.b Running DIA-NN to perform LFQ
example: 
```
#!/bin/bash 
set -ex 
mkdir -p temp-DIANN 
mkdir -p out-DIANN 
nice -19 PATH/TO/diann-linux \
--fasta-search --fasta DATABASE.fasta \
--f 20251016_002_S1030253_MACE_snow_test.mzML  \
--threads 64 --qvalue 0.01 --matrices --predictor \
--met-excision --cut K*,R* --min-pep-len 6 --max-pep-len 30 --smart-profiling \
--var-mods 1 --var-mod UniMod:35,15.994915,M \
--min-pr-charge 2 --max-pr-charge 3 \
--min-pr-mz 400 --max-pr-mz 1500 \
--verbose 1 --unimod4 \
--missed-cleavages 1 \
--mass-acc 15 \
--mass-acc-ms1 10 \
--reanalyse --pg-level 0 \
--out-lib out-DIANN/REPORT_LIB.tsv --out-lib-copy \
--temp temp-DIANN \
--out out-DIANN/REPORT_NAME.tsv \
| tee diann.log.txt
```

# 3. Extract protein intensities with PickedGroupFDR
For Sage: 
```
python3 -u -m picked_group_fdr \
   --protein_group_fdr_threshold 0.05 \
   --fasta DATABASE \
   --sage_results results.sage.tsv \
   --sage_lfq_tsv lfq.tsv \
   --protein_groups_out combined_protein.tsv \
   --output_format fragpipe \
   --do_quant \
   --lfq_min_peptide_ratios 1 \
   --methods sage
```

For DIA-NN:
```
python3 -u -m picked_group_fdr \
   --fasta DATABASE \
   --fasta_use_uniprot_id \
   --diann_reports ./diann/report.parquet \
   --protein_groups_out combined_protein.tsv \
   --output_format fragpipe \
   --do_quant \
   --lfq_min_peptide_ratios 1 \
   --methods diann
```
