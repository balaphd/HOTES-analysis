
########################## BACTERIA ##########################
##uses the split_libraries_fastq command to convert from fastq to fasta with quality filtering while definding sample names for each sample: resequenced individuals were given the suffix “R”. Because the command was done per sample, a dummy mapping file was used. ##
##combine the sequences##
cat bact_slout_q20/slout_single_sample_q20_16s*/seqs.fna > final_demultiplexed_seqs_16S.fasta 
##identify (using usearch wrapped within QIIME) then exlucde chimeras##
identify_chimeric_seqs.py -i final_demultiplexed_seqs_16S.fasta -m usearch61 -o bact_usearch_checked_chimeras/ -r greengenes/gg_13_8_otus/rep_set/97_otus.fasta
filter_fasta.py -f final_demultiplexed_seqs_16S.fasta -o final_chimeric_rmv_seqs_16S.fasta -s bact_usearch_checked_chimeras/chimeras.txt -n
##cluster sequences into OTUs, first by closed-reference against the greengenes database then by de novo clustering at 100% subsampling ##
pick_open_reference_otus.py -i final_chimeric_rmv_seqs_16S.fasta  -r greengenes/gg_13_8_otus/rep_set/97_otus.fasta -o bact_open_ref_picked_otus -a -O 4 --percent_subsample 1
## filter spurious OTUs from total sequences
filter_otus_from_otu_table.py -i bact_open_ref_picked_otus/otu_table_mc2_w_tax_no_pynast_failures.biom -o bact_otu_table_mc2_w_tax_no_pynast_failures_min00001.biom --min_count_fraction 0.00001
## looking at unassigned sequences (n=17), they mostly BLAST to other uncultured bacterium; given their low number and their likelihood of being undescribed bacteria, they were left in. 
## remove all chloroplasts/mitochondria then all empty OTUs
filter_taxa_from_otu_table.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001.biom -o bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -n  c__Chloroplast,f__mitochondria
(confirmed those sequences were removed)
(grep -o 'Chloroplast' bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom | wc -l)
       0
(grep -o 'mitochondria' bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom | wc -l)
       0
filter_otus_from_otu_table.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -o bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro.biom -n 1
## rarefaction value: was selected based on the resonable lowest cutoff 
## remove duplicates, controls, mock, etc
filter_samples_from_otu_table.py -i bact_tab_after_filtering_no_unassigned_no_chloro.biom -o bact_tab_after_filtering_no_unassigned_no_chloro_final_samples.biom --sample_id_fp bact_final_samples.txt
filter_otus_from_otu_table.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -o bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -n 1
biom convert -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -o bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.txt --header-key taxonomy -b
## final bacteria OTU table 
bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom
### make tree of final OTUs
filter_fasta.py -f bact_open_ref_picked_otus/rep_set.fna -o bact_final_OTUs.fna -b bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom
parallel_align_seqs_pynast.py -i bact_final_OTUs.fna -o bact_final_OTUs_pynast_aligned_seqs/ -O 6
filter_alignment.py -o bact_final_OTUs_pynast_aligned_seqs/ -i bact_final_OTUs_pynast_aligned_seqs/bact_final_OTUs_aligned.fasta --allowed_gap_frac 0.999999 --threshold 3.0 --suppress_lane_mask_filter
make_phylogeny.py -i bact_final_OTUs_pynast_aligned_seqs/bact_final_OTUs_aligned_pfiltered.fasta -o bact_final_OTUs.tre --root_method tree_method_default --tree_method fasttree
########################## FUNGI ###########################


cat fung_slout_q20/slout_single_sample_q20_ITS*/seqs.fna > final_demultiplexed_seqs_ITS_all.fasta
identify_chimeric_seqs.py -i final_demultiplexed_seqs_ITS.fasta -m usearch61 -o fung_usearch_checked_chimeras/ -r reference_db/its_12_11_otus/rep_set/97_otus.fasta
filter_fasta.py -f final_demultiplexed_seqs_ITS_all.fasta -o final_chimeric_rmv_seqs_ITS_step1.fasta -s fung_usearch_checked_chimeras/chimeras.txt -n
pick_open_reference_otus.py -i final_chimeric_rmv_seqs_ITS.fasta -r reference_db/its_12_11_otus/rep_set/97_otus.fasta -o fungi_open_ref_picked_otus --prefilter_percent_id 0.0 --suppress_align_and_tree -a -O 6 --percent_subsample 1
## tried several ways of assigning taxonomy 
assign_taxonomy.py -i fungi_open_ref_picked_otus/rep_set.fna -t reference_db/sh_qiime_release_02.03.2015/sh_taxonomy_qiime_ver7_97_02.03.2015.txt -r reference_db/sh_qiime_release_02.03.2015/sh_refs_qiime_ver7_97_02.03.2015.fasta -o fung_qiime_assigned_tax_sortmerna -m sortmerna
## filter spurious OTUs from total sequences
filter_otus_from_otu_table.py -i fungi_open_ref_picked_otus/otu_table_mc2.biom -o fung_otu_table_mc2_min00001.biom --min_count_fraction 0.00001
##biom taxonomy to the biom file 
biom add-metadata -i fung_otu_table_mc2_min00001.biom -o fung_otu_table_mc2_min00001_w_sortmerna_tax.biom --observation-metadata-fp fung_qiime_assigned_tax_sortmerna/rep_set_tax_assignments_header.txt --sc-separated taxonomy
biom convert -i fung_otu_table_mc2_min00001_w_sortmerna_tax.biom -o fung_otu_table_mc2_min00001_w_sortmerna_tax.txt --header-key taxonomy -b
## remove duplicates, controls, mock, samples with too low read count 
filter_samples_from_otu_table.py -i fung_otu_table_mc2_min00001_w_sortmerna_tax.biom -o fung_otu_table_mc2_min00001_w_sortmerna_tax_final_samples.biom --sample_id_fp fung_final_samples.txt
biom convert -i fung_otu_table_mc2_min00001_w_sortmerna_tax_final_samples.biom -o fung_otu_table_mc2_min00001_w_sortmerna_tax_final_samples.txt --header-key taxonomy -b
biom summarize-table -i fung_otu_table_mc2_min00001_w_sortmerna_tax_final_samples.biom -o fung_otu_table_mc2_min00001_w_sortmerna_tax_final_samples_summary.txt
## final fungal OTU table 
fung_otu_table_mc2_min00001_w_rdp_tax_final_samples.biom 





########## community analyses###############

## alpha richness across samples ##
alpha_rarefaction.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -m map.hotes.txt -o alpha_rare_e1000_hotes -a -t rep_set.tre -e 1000 -O 6 -p alpha_params.txt
## beta diversity across samples ##
beta_diversity_through_plots.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -m map.hotes.txt -o bdiv_e1000_hotes -a -t rep_set.tre -e 1000 -O 6 -p beta_params.txt
#############################################
summarize_taxa.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -o bact_tax/
summarize_taxa.py -i bact_otu_table_mc2_w_tax_no_pynast_failures_min00001_nochloro_final_samples.biom -o bact_tax_abs/ -a

####################################################
