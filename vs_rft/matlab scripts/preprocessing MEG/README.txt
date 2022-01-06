Script order

1 trial structure
	a_merge_edf_fif_mat.mat: identifies and aligns trials found in MEG, responses (mat file) and edf file
	resulting trial structure to be used for following analyses

2 soi pipeline: identifies soi with sign tagging response
	a_ICA: ICA using trial structure identified above (function)
	a2_ICA_identify_comp: identification of bad components by user (script)
	b_find_soi_stats: loads in data as per structure identified in 1), 
			  suppresses bad components identifed in a2_[...] in 2 a2), 
			  runs cluster-perm statistical test on stim vs baseline,
			  returns significant sensors
