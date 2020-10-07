From Frank's shuffled Z-score matrix, the process is:
1. Remove tissues with <9 cell lines. Leaves 550 cell lines.
2. Calculate gene-by-gene correlation matrix of all genes, 550 cell lines (cc_zscore_matrix)
3. Calculate correlation matrix of AML-only, 15 cell lines (cc_matrix_aml15cells)
4. Calculate correlation matrix of 530 cell lines, excluding blood lineage (cc_matrix_exBlood)
5. Calculate difference in correlation for each pair (delta-cc-list)
6. Extract information from data!


Differences in correlations are not always meaningful. The get_sig_of_corrs script(s) attempt
to empirically measure the significance of delta-PCC. Since there are 15 AML cell lines, the
script randomly selects 15 cell lines and calculates all gene-gene correlations, keeping track
of whether the random correlation is greater (or less) than the AML-only PCC. This process is
repeated 1000x and an empirical P-value up to 0.001 is determined.