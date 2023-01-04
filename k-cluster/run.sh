kclust -mode rmsd -radius 9 -centroid Abeta40_filelist > Abeta40_cluster.dat
perl analyse_cluster_percentage.pl Abeta40_cluster.dat ./re_kclust/Abeta40_processed.dat 90.0

