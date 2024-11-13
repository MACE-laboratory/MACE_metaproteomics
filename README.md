# Metaproteomics pipeline for the MACE laboratory
Massimo Bourquin, 2024

# 1. Raw to mzML conversion
example: 
```
thermorawfileparser -i SAMPLE.raw  -N -f 1 -g
```

# 2. Running Sage to perform LFQ
example: 
```
sage config.json --batch-size 12
```

# 3. Extract protein intensities with PickedGroupFDR
example: 
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
