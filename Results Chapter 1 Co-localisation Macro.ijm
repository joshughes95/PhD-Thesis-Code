name = getTitle();
selectWindow(name);
run("8-bit");
run("Subtract Background...", "rolling=50 stack");
run("Split Channels");

run("Coloc 2", "channel_1=[C1-"+name+"] channel_2=[C2-"+name+"] roi_or_mask=<None> threshold_regression=Costes show_save_pdf_dialog li_histogram_channel_1 li_histogram_channel_2 li_icq spearman's_rank_correlation manders'_correlation kendall's_tau_rank_correlation 2d_intensity_histogram costes'_significance_test psf=3 costes_randomisations=10");
