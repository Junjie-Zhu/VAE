#!/usr/bin/perl
for($j=0;$j<500;$j+=1){
	$i = $j * 49 + 1;
	print	join($i, "./pdbs_to_use/Abeta40/Abeta40_", ".pdb\n");
}
