die "Usage: perl $0 <ClusterFile> <OUTPUT_Prefix> <Percentage>\n" unless @ARGV == 3;
$Input = $ARGV[0];
$Output = $ARGV[1];
$MaxC = $ARGV[2];

open (IN,$Input) || die "File $Input open failed!!\n";
$totalfile = 0;
$i = -1;
while (<IN>){
	if (/^\#Cluster \d+/){
		$i ++;
		$lowest[$i] = 2000;
		while (1 != 0){
			chomp ($tmp = <IN>);
			last if $tmp =~ /^\#/;
			$ClusterSize[$i] ++;
			$totalfile ++;
			@temp = split /\s+/,$tmp;
			if ($temp[2] < $lowest[$i]){
				$lowest[$i] = $temp[2];
				$Pick[$i] = $temp[1];
			}
		}
	}
}
close IN;
$ClusterNum = $i+1;
print "$ClusterNum clusters found!\nSorting...\n";

foreach $i(0..$ClusterNum -2){
	foreach $j($i+1 .. $ClusterNum-1){
		if ($ClusterSize[$i] < $ClusterSize[$j]){
			($ClusterSize[$i],$ClusterSize[$j]) = ($ClusterSize[$j],$ClusterSize[$i]);
			($lowest[$i],$lowest[$j]) = ($lowest[$j],$lowest[$i]);
			($Pick[$i], $Pick[$j]) = ($Pick[$j], $Pick[$i]);
		}
	}
}

open OUT,"> $Output\_Summary";
print OUT "Total Structures: $totalfile\n";
print OUT "Total Clusters: $ClusterNum\n\n";
$TTOCP = 0;
$outputNum = 0;
while ($TTOCP < $MaxC){
	$j = $outputNum +1;
	printf OUT "%d\t%.3f\%\n",$outputNum+1, $ClusterSize[$outputNum]/$totalfile*100;
	$TTOCP += $ClusterSize[$outputNum]/$totalfile*100;
	print $outputNum,"\t$ClusterSize[$outputNum]\t$lowest[$outputNum]\t$Pick[$outputNum]\n";
	`cp $Pick[$outputNum] $Output\_$j.pdb`;
	$outputNum ++;
}

printf OUT "TOP %d Clusters Total: %.3f\%\n", $outputNum, $TTOCP;
close OUT;
open OUT,"> $Output\_Alignment.pml";
foreach $i(1..$outputNum){
	print OUT "load $Output\_$i.pdb\n";
}
foreach $i(2..$outputNum){
	print OUT "align $Output\_$i\/\/\/\/CA, $Output\_1\n";
}
close OUT;
