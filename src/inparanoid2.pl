#! /usr/bin/perl

###RAM UPDATES
#inparanoid2.pl - checks for existence of BLAST files before running time-consuming BLAST.
#This is important for "inside" scores, instead of running that process "n" times,
#they are only run once...
#Also massive amount of stripping out HTML and bootstrap functionality
#removing do_blast and just generating diamond files with DiamondParser.py


use File::Copy;
my $usage =" Usage: inparanoid.pl <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <path to outdir>";

$run_inparanoid = 1;
$table = 1;       # Print tab-delimited table of orthologs to file "table.txt" #

# Algorithm parameters:
# Default values should work without problems.
# MAKE SURE, however, that the score cutoff here matches what you used for BLAST!
$score_cutoff = 40;    # In bits. Any match below this is ignored             #
$outgroup_cutoff = 50; # In bits. Outgroup sequence hit must be this many bits#
                       # stronger to reject best-best hit between A and B     #
$conf_cutoff = 0.05;   # Include in-paralogs with this confidence or better   #
$group_overlap_cutoff = 0.5; # Merge groups if ortholog in one group has more #
                             # than this confidence in other group            #
$grey_zone = 0;  # This many bits signifies the difference between 2 scores   #
$debug = 0;      # Print debugging messages or not. Levels 0,1,2 and 4 exist  #

my $seq_overlap_cutoff = 0.5; 		# Match area should cover at least this much of longer sequence. Match area is defined as area from start of
					# first segment to end of last segment, i.e segments 1-10 and 90-100 gives a match length of 100.
my $segment_coverage_cutoff = 0.25; 	# Actually matching segments must cover this much of longer sequence.
					# For example, segments 1-10 and 90-100 gives a total length of 20.

###############################################################################
# No changes should be required below this line                               #
###############################################################################

if (!@ARGV){
    print STDERR $usage;
    exit 1;
}

if ((@ARGV < 2) and ($run_inparanoid)){
    print STDERR "\n When \$run_inparanoid=1, at least two distinct FASTA files have to be specified.\n";

    print STDERR $usage;
    exit 1;
}

if ((!$run_blast) and (!$run_inparanoid)){
    print STDERR "run_blast or run_inparanoid has to be set!\n";
    exit 1;
}

# Input files:
$fasta_seq_fileA = "$ARGV[2]" . "faa/" . "$ARGV[0]" . ".faa";
$fasta_seq_fileB = "$ARGV[2]" . "faa/" . "$ARGV[1]" . ".faa";

my $blast_outputAB = "$ARGV[2]" . "out/" . "$ARGV[0]" . "." . "$ARGV[1]" . ".out";
my $blast_outputBA = "$ARGV[2]" . "out/" . "$ARGV[1]" . "." . "$ARGV[0]" . ".out";
my $blast_outputAA = "$ARGV[2]" . "out/" . "$ARGV[0]" . "." . "$ARGV[0]" . ".out";
my $blast_outputBB = "$ARGV[2]" . "out/" . "$ARGV[1]" . "." . "$ARGV[1]" . ".out";


my %idA;        # Name -> ID combinations for species 1
my %idB;        # Name -> ID combinations for species 2
my @nameA;      # ID -> Name combinations for species 1
my @nameB;      # ID -> Name combinations for species 2
my %scoreAB;    # Hashes with pairwise BLAST scores (in bits)
my %scoreBA;
my %scoreAA;
my %scoreBB;
my @hitnAB;     # 1-D arrays that keep the number of pairwise hits
my @hitnBA;
my @hitnAA;
my @hitnBB;
my @hitAB;      # 2-D arrays that keep the actual matching IDs
my @hitBA;
my @hitAA;
my @hitBB;
my @besthitAB;  # IDs of best hits in other species (may contain more than one ID)
my @besthitBA;  # IDs of best hits in other species (may contain more than one ID)
my @bestscoreAB; # best match A -> B
my @bestscoreBA; # best match B -> A
my @ortoA;      # IDs of ortholog candidates from species A
my @ortoB;      # IDs of ortholog candidates from species B
my @ortoS;      # Scores between ortoA and ortoB pairs
my @paralogsA;  # List of paralog IDs in given cluster
my @paralogsB;  # List of paralog IDs in given cluster
my @confPA;     # Confidence values for A paralogs
my @confPB;     # Confidence values for B paralogs
my @confA;      # Confidence values for orthologous groups
my @confB;      # Confidence values for orthologous groups
my $prev_time = 0;

#################################################
# Assign ID numbers for species A
#################################################
open A, "$fasta_seq_fileA" or die "File A with sequences in FASTA format is missing
Usage $0 <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <FASTAFILE with sequences of species C>\n";
$id = 0;
while (<A>){
    if(/^\>/){
	++$id;
	chomp;
	s/\>//;
	@tmp = split(/\s+/);
	#$name = substr($tmp[0],0,25);
	$name = $tmp[0];
	$idA{$name} = int($id);
	$nameA[$id] = $name;
    }
}
close A;
$A = $id;

if (@ARGV >= 2) {
#################################################
# Assign ID numbers for species B
#################################################
    open B, "$fasta_seq_fileB" or die "File B with sequences in FASTA format is missing
Usage $0 <FASTAFILE with sequences of species A> <FASTAFILE with sequences of species B> <FASTAFILE with sequences of species C>\n";
    $id = 0;
    while (<B>){
	if(/^\>/){
	    ++$id;
	    chomp;
	    s/\>//;
	    @tmp = split(/\s+/);
	    #$name = substr($tmp[0],0,25);
	    $name = $tmp[0];
	    $idB{$name} = int($id);
	    $nameB[$id] = $name;
	}
    }
    $B = $id;
    close B;

}

if ($run_inparanoid){
#################################################
# Read in best hits from blast output file AB
#################################################
    $count = 0;
    open AB, "$blast_outputAB" or die "Blast output file A->B is missing\n";
    $old_idQ = 0;
    while (<AB>){
	chomp;
	@Fld = split(/\s+/);    # Get query, match and score

	if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
		print STDERR "AB ok\n";
	    }
	    next;
	}

	$q = $Fld[0];
	$m = $Fld[1];
	$idQ = $idA{$q}; # ID of query sequence
	$idM = $idB{$m}; # ID of match sequence
	$score = $Fld[2];

	next if (!overlap_test(@Fld));

	# Score must be equal to or above cut-off
	next if ($score < $score_cutoff);

	if(!$count || $q ne $oldq){
	    print "Match $m, score $score, ID for $q is missing\n" if ($debug == 2 and !(exists($idA{$q})));
	    $hitnAB[$idA{$oldq}] = $hit if($count); # Record number of hits for previous query
	    $hit = 0;
	    ++$count;
	    $oldq = $q;
	}
	++$hit;
	$hitAB[$idQ][$hit] = int($idM);
#	printf ("hitAB[%d][%d] = %d\n",$idQ,$hit,$idM);
	$scoreAB{"$idQ:$idM"} = $score;
	$scoreBA{"$idM:$idQ"} = $score_cutoff; # Initialize mutual hit score - sometimes this is below score_cutoff
	$old_idQ = $idQ;
#    }
    }
    $hitnAB[$idQ] = $hit; # For the last query
#printf ("hitnAB[1] = %d\n",$hitnAB[1]);
#printf ("hitnAB[%d] = %d\n",$idQ,$hit);
    close AB;

#################################################
# Read in best hits from blast output file BA
#################################################
    $count = 0;
    open BA, "$blast_outputBA" or die "Blast output file B->A is missing\n";
    $old_idQ = 0;
    while (<BA>){
	chomp;
	@Fld = split(/\s+/);    # Get query, match and score

	if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
		print STDERR "BA ok\n";
	    }
	    next;
	}

	$q = $Fld[0];
	$m = $Fld[1];
	$idQ = $idB{$q};
	$idM = $idA{$m};
	$score = $Fld[2];

	next if (!overlap_test(@Fld));

	next if ($score < $score_cutoff);

	if(!$count || $q ne $oldq){
	    print "ID for $q is missing\n" if ($debug == 2 and (!exists($idB{$q})));
	    $hitnBA[$idB{$oldq}] = $hit if($count); # Record number of hits for previous query
	    $hit = 0;
	    ++$count;
	    $oldq = $q;
	}
	++$hit;
	$hitBA[$idQ][$hit] = int($idM);
	$scoreBA{"$idQ:$idM"} = $score;
	$scoreAB{"$idM:$idQ"} = $score_cutoff if ($scoreAB{"$idM:$idQ"} < $score_cutoff); # Initialize missing scores
	$old_idQ = $idQ;
    }
    $hitnBA[$idQ] = $hit; # For the last query
    close BA;
##################### Equalize AB scores and BA scores ##########################

    	foreach my $key (keys %scoreAB) {

		my ($a, $b) = split(':', $key);
		my $key2 = $b . ':' . $a;

		# If debugg mod is 5 and the scores A-B and B-A are unequal
	   	 # the names of the two sequences and their scores are printed
	    	if ($scoreAB{$key} != $scoreBA{$key2}){
			printf ("%-20s\t%-20s\t%d\t%d\n",$nameA[$a], $nameB[$b], $scoreAB{$key}, $scoreBA{$key2}) if ($debug == 5);
	    	}

		# Set score AB and score BA to the mean of scores AB and BA.
	  	# The final score is saved as an integer so .5 needs to be added to avoid rounding errors
	    	$scoreAB{$key} = $scoreBA{$key2} = int(($scoreAB{$key} + $scoreBA{$key2})/2.0 +.5);
	}

##################### Re-sort hits, besthits and bestscore #######################
    for $idA(1..$A){

	next if (!($hitnAB[$idA]));
	for $hit (1..($hitnAB[$idA]-1)){ # Sort hits by score
	    while($scoreAB{"$idA:$hitAB[$idA][$hit]"} < $scoreAB{"$idA:$hitAB[$idA][$hit+1]"}){
		$tmp = $hitAB[$idA][$hit];
		$hitAB[$idA][$hit] = $hitAB[$idA][$hit+1];
		$hitAB[$idA][$hit+1] = $tmp;
		--$hit if ($hit > 1);
	    }
	}
	$bestscore = $bestscoreAB[$idA] = $scoreAB{"$idA:$hitAB[$idA][1]"};
	$besthitAB[$idA] = $hitAB[$idA][1];
	for $hit (2..$hitnAB[$idA]){
	    if ($bestscore - $scoreAB{"$idA:$hitAB[$idA][$hit]"} <= $grey_zone){
		$besthitAB[$idA] .= " $hitAB[$idA][$hit]";
	    }
	    else {
		last;
	    }
	}
	undef $is_besthitAB[$idA]; # Create index that we can check later
	grep (vec($is_besthitAB[$idA],$_,1) = 1, split(/ /,$besthitAB[$idA]));


    }

    for $idB(1..$B){
#    print "Loop index = $idB\n";
	next if (!($hitnBA[$idB]));
	for $hit (1..($hitnBA[$idB]-1)){
# Sort hits by score
	    while($scoreBA{"$idB:$hitBA[$idB][$hit]"} < $scoreBA{"$idB:$hitBA[$idB][$hit+1]"}){
		$tmp = $hitBA[$idB][$hit];
		$hitBA[$idB][$hit] = $hitBA[$idB][$hit+1];
		$hitBA[$idB][$hit+1] = $tmp;
		--$hit if ($hit > 1);
	    }
	}
	$bestscore = $bestscoreBA[$idB] = $scoreBA{"$idB:$hitBA[$idB][1]"};
	$besthitBA[$idB] = $hitBA[$idB][1];
	for $hit (2..$hitnBA[$idB]){
	    if ($bestscore - $scoreBA{"$idB:$hitBA[$idB][$hit]"} <= $grey_zone){
		$besthitBA[$idB] .= " $hitBA[$idB][$hit]";
	    }
	    else {last;}
	}
	undef $is_besthitBA[$idB]; # Create index that we can check later
	grep (vec($is_besthitBA[$idB],$_,1) = 1, split(/ /,$besthitBA[$idB]));
#    printf ("besthitBA[%d] = %d\n",$idA,$besthitAB[$idA]);
    }


######################################################
# Now find orthologs:
######################################################
    $o = 0;

    for $i(1..$A){   # For each ID in file A
	if (defined $besthitAB[$i]){
	    @besthits = split(/ /,$besthitAB[$i]);
	    for $hit (@besthits){
		if (vec($is_besthitBA[$hit],$i,1)){
		    ++$o;
		    $ortoA[$o] = $i;
		    $ortoB[$o] = $hit;
		    $ortoS[$o] = $scoreAB{"$i:$hit"}; # Should be equal both ways
#	    --$o if ($ortoS[$o] == $score_cutoff); # Ignore orthologs that are exactly at score_cutoff
		    print "Accept! " if ($debug == 2);
		}
		else {print "        " if ($debug == 2);}
		printf ("%-20s\t%d\t%-20s\t", $nameA[$i], $bestscoreAB[$i], $nameB[$hit]) if ($debug == 2);
		print "$bestscoreBA[$hit]\t$besthitBA[$hit]\n" if ($debug == 2);
	    }
	}
    }
    print "$o ortholog candidates detected\n" if ($debug);
#####################################################
# Sort orthologs by ID and then by score:
#####################################################

    # Create an array used to store the position each element shall have in the final array
    # The elements are initialized with the position numbers
    my @position_index_array = (1..$o);

    # Sort the position list according to id
    my @id_sorted_position_list = sort { ($ortoA[$a]+$ortoB[$a]) <=> ($ortoA[$b] + $ortoB[$b]) } @position_index_array;

    # Sort the list according to score
    my @score_id_sorted_position_list = sort { $ortoS[$b] <=> $ortoS[$a] } @id_sorted_position_list;

    # Create new arrays for the sorted information
    my @new_ortoA;
    my @new_ortoB;
    my @new_orthoS;

   # Add the information to the new arrays in the orer specifeid by the index array
   for (my $index_in_list = 0; $index_in_list < scalar @score_id_sorted_position_list; $index_in_list++) {


	my $old_index = $score_id_sorted_position_list[$index_in_list];
	$new_ortoA[$index_in_list + 1] = $ortoA[$old_index];
	$new_ortoB[$index_in_list + 1] = $ortoB[$old_index];
	$new_ortoS[$index_in_list + 1] = $ortoS[$old_index];
   }

    @ortoA = @new_ortoA;
    @ortoB = @new_ortoB;
    @ortoS = @new_ortoS;

    @all_ortologsA = ();
    @all_ortologsB = ();
    for $i(1..$o){
	push(@all_ortologsA,$ortoA[$i]); # List of all orthologs
	push(@all_ortologsB,$ortoB[$i]);
    }
    undef $is_ortologA; # Create index that we can check later
    undef $is_ortologB;
    grep (vec($is_ortologA,$_,1) = 1, @all_ortologsA);
    grep (vec($is_ortologB,$_,1) = 1, @all_ortologsB);


################################
# Read inside scores from AA
################################
    $count = 0;
    $max_hit = 0;
    open AA, "$blast_outputAA" or die "Blast output file A->A is missing\n";
    while (<AA>) {
	chomp;                  # strip newline

	@Fld = split(/\s+/);    # Get query and match names

	if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
		print STDERR "AA ok\n";
	    }
	    next;
	}

	$q = $Fld[0];
	$m = $Fld[1];
	$score = $Fld[2];
	next unless (vec($is_ortologA,$idA{$q},1));

	next if (!overlap_test(@Fld));

	next if ($score < $score_cutoff);

	if(!$count || $q ne $oldq){ # New query
	    $max_hit = $hit if ($hit > $max_hit);
	    $hit = 0;
	    $oldq = $q;
	}
	++$hit;
	++$count;
	$scoreAA{"$idA{$q}:$idA{$m}"}  = int($score + 0.5);
	$hitAA[$idA{$q}][$hit] = int($idA{$m});
	$hitnAA[$idA{$q}] = $hit;
    }
    close AA;

################################
# Read inside scores from BB
################################
    $count = 0;
    open BB, "$blast_outputBB" or die "Blast output file B->B is missing\n";
    while (<BB>) {
	chomp;                  # strip newline

	@Fld = split(/\s+/);    # Get query and match names

	if( scalar @Fld < 9 ){
	    if($Fld[0]=~/done/){
		print STDERR "BB ok\n";
	    }
	    next;
	}

	$q = $Fld[0];
	$m = $Fld[1];
	$score = $Fld[2];
	next unless (vec($is_ortologB,$idB{$q},1));

	next if (!overlap_test(@Fld));

	next if ($score < $score_cutoff);

	if(!$count || $q ne $oldq){ # New query
	    $max_hit = $hit if ($hit > $max_hit);
	    $oldq = $q;
	    $hit = 0;
	}
	++$count;
	++$hit;
	$scoreBB{"$idB{$q}:$idB{$m}"} = int($score + 0.5);
	$hitBB[$idB{$q}][$hit] = int($idB{$m});
	$hitnBB[$idB{$q}] = $hit;
    }
    close BB;

    print "Maximum number of hits per sequence was $max_hit\n" if ($debug);
#####################################################
# Find paralogs:
#####################################################
    for $i(1..$o){
	$merge[$i] = 0;
	next if($del[$i]); # If outgroup species was closer to one of the seed orthologs
	$idA = $ortoA[$i];
	$idB = $ortoB[$i];
	local @membersA = ();
	local @membersB = ();

	undef $is_paralogA[$i];
	undef $is_paralogB[$i];

	print "$i: Ortholog pair $nameA[$idA] and $nameB[$idB]. $hitnAA[$idA] hits for A and $hitnBB[$idB] hits for B\n"  if ($debug);
	# Check if current ortholog is already clustered:
	for $j(1..($i-1)){
	    # Overlap type 1: Both orthologs already clustered here -> merge
	    if ((vec($is_paralogA[$j],$idA,1)) and (vec($is_paralogB[$j],$idB,1))){
		$merge[$i] = $j;
		print "Merge CASE 1: group $i ($nameB[$idB]-$nameA[$idA]) and $j ($nameB[$ortoB[$j]]-$nameA[$ortoA[$j]])\n" if ($debug);
		last;
	    }
	    # Overlap type 2: 2 competing ortholog pairs -> merge
	    elsif (($ortoS[$j] - $ortoS[$i] <= $grey_zone)
		   and (($ortoA[$j] == $ortoA[$i]) or ($ortoB[$j] == $ortoB[$i]))
#       and ($paralogsA[$j])
		   ){ # The last condition is false if the previous cluster has been already deleted
		$merge[$i] = $j;
		print "Merge CASE 2: group $i ($nameA[$ortoA[$i]]-$nameB[$ortoB[$i]]) and $j ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
		last;
	    }
	    # Overlap type 3: DELETE One of the orthologs belongs to some much stronger cluster -> delete
	    elsif (((vec($is_paralogA[$j],$idA,1)) or (vec($is_paralogB[$j],$idB,1))) and
		   ($ortoS[$j] - $ortoS[$i] > $score_cutoff)){
		print "Delete CASE 3: Cluster $i -> $j, score $ortoS[$i] -> $ortoS[$j], ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
		$merge[$i]= -1; # Means - do not add sequences to this cluster
		$paralogsA[$i] = "";
		$paralogsB[$i] = "";
		last;
	    }
	    # Overlap type 4: One of the orthologs is close to the center of other cluster
	    elsif (((vec($is_paralogA[$j],$idA,1)) and ($confPA[$idA] > $group_overlap_cutoff)) or
		   ((vec($is_paralogB[$j],$idB,1)) and ($confPB[$idB] > $group_overlap_cutoff))){
		print "Merge CASE 4: Cluster $i -> $j, score $ortoS[$i] -> $ortoS[$j], ($nameA[$ortoA[$j]]-$nameB[$ortoB[$j]])\n" if ($debug);
		$merge[$i] = $j;
		last;
	    }
	    # Overlap type 5:
	    # All clusters that were overlapping, but not catched by previous "if" statements will be DIVIDED!
	}
	next if ($merge[$i] < 0); # This cluster should be deleted
##### Check for paralogs in A
	$N = $hitnAA[$idA];
	for $j(1..$N){
	    $hitID = $hitAA[$idA][$j]; # hit of idA
#      print "Working with $nameA[$hitID]\n" if ($debug == 2);
	    # Decide whether this hit is inside the paralog circle:
	    if ( ($idA == $hitID) or ($scoreAA{"$idA:$hitID"} >= $bestscoreAB[$idA]) and
		($scoreAA{"$idA:$hitID"} >= $bestscoreAB[$hitID])){
		if ($debug == 2){
		    print "   Paralog candidates: ";
		    printf ("%-20s: %-20s", $nameA[$idA], $nameA[$hitID]);
		    print "\t$scoreAA{\"$idA:$hitID\"} : $bestscoreAB[$idA] : $bestscoreAB[$hitID]\n";
		}
		$paralogs = 1;
		if ($scoreAA{"$idA:$idA"} == $ortoS[$i]){
		    if ($scoreAA{"$idA:$hitID"} == $scoreAA{"$idA:$idA"}){
			$conf_here = 1.0; # In the center
		    }
		    else{
			$conf_here = 0.0; # On the border
		    }
		}
		else {
		    $conf_here = ($scoreAA{"$idA:$hitID"} - $ortoS[$i]) /
			($scoreAA{"$idA:$idA"} - $ortoS[$i]);
		}
		# Check if this paralog candidate is already clustered in other clusters
		for $k(1..($i-1)){
		    if (vec($is_paralogA[$k],$hitID,1)){ # Yes, found in cluster $k
			if($debug == 2){
			    print "      $nameA[$hitID] is already in cluster $k, together with:";
			    print " $nameA[$ortoA[$k]] and $nameB[$ortoB[$k]] ";
			    print "($scoreAA{\"$ortoA[$k]:$hitID\"})";
			}
			if (($confPA[$hitID] >= $conf_here) and
			    ($j != 1)){ # The seed ortholog CAN NOT remain there
			    print " and remains there.\n" if ($debug == 2);
			    $paralogs = 0; # No action
			}
			else { # Ortholog of THIS cluster is closer than ortholog of competing cluster $k
			    print " and should be moved here!\n" if ($debug == 2); # Remove from other cluster, add to this cluster
			    @membersAK = split(/ /, $paralogsA[$k]); # This array contains IDs
			    $paralogsA[$k] = "";# Remove all paralogs from cluster $k
				@tmp = ();
			    for $m(@membersAK){
				push(@tmp,$m) if ($m != $hitID); # Put other members back
			    }
			    $paralogsA[$k] = join(' ',@tmp);
			    undef $is_paralogA[$k]; # Create index that we can check later
			    grep (vec($is_paralogA[$k],$_,1) = 1, @tmp);
			}
			last;
		    }
		}
		next if (! $paralogs); # Skip that paralog - it is already in cluster $k
		push (@membersA,$hitID); # Add this hit to paralogs of A
	    }
	}
	# Calculate confidence values now:
	@tmp = ();
	for $idP (@membersA){ # For each paralog calculate conf value
	    if($scoreAA{"$idA:$idA"} == $ortoS[$i]){
		if ($scoreAA{"$idA:$idP"} == $scoreAA{"$idA:$idA"}){
		    $confPA[$idP] = 1.00;
		}
		else{
		    $confPA[$idP] = 0.00;
		}
	    }
	    else{
		$confPA[$idP] = ($scoreAA{"$idA:$idP"} - $ortoS[$i]) /
		    ($scoreAA{"$idA:$idA"} - $ortoS[$i]);
	    }
	    push (@tmp, $idP) if ($confPA[$idP] >= $conf_cutoff); # If one wishes to use only significant paralogs
	}
	@membersA = @tmp;
	########### Merge if necessary:
	if ($merge[$i] > 0){ # Merge existing cluster with overlapping cluster
	    @tmp = split(/ /,$paralogsA[$merge[$i]]);
	    for $m (@membersA){
		push (@tmp, $m)  unless (vec($is_paralogA[$merge[$i]],$m,1));
	    }
	    $paralogsA[$merge[$i]] = join(' ',@tmp);
	    undef $is_paralogA[$merge[$i]];
	    grep (vec($is_paralogA[$merge[$i]],$_,1) = 1, @tmp); # Refresh index of paralog array
	}
	######### Typical new cluster:
	else{  # Create a new cluster
	    $paralogsA[$i] = join(' ',@membersA);
	    undef $is_paralogA; # Create index that we can check later
	    grep (vec($is_paralogA[$i],$_,1) = 1, @membersA);
	}
##### The same procedure for species B:
	$N = $hitnBB[$idB];
	for $j(1..$N){
	    $hitID = $hitBB[$idB][$j];
#      print "Working with $nameB[$hitID]\n" if ($debug == 2);
	    if ( ($idB == $hitID) or ($scoreBB{"$idB:$hitID"} >= $bestscoreBA[$idB]) and
		($scoreBB{"$idB:$hitID"} >= $bestscoreBA[$hitID])){
		if ($debug == 2){
		    print "   Paralog candidates: ";
		    printf ("%-20s: %-20s", $nameB[$idB], $nameB[$hitID]);
		    print "\t$scoreBB{\"$idB:$hitID\"} : ";
		    print "$bestscoreBA[$idB] : $bestscoreBA[$hitID]\n";
		}
		$paralogs = 1;
		if ($scoreBB{"$idB:$idB"} == $ortoS[$i]){
		    if ($scoreBB{"$idB:$hitID"} == $scoreBB{"$idB:$idB"}){
			$conf_here = 1.0;
		    }
		    else{
			$conf_here = 0.0;
		    }
		}
		else{
		    $conf_here = ($scoreBB{"$idB:$hitID"} - $ortoS[$i]) /
			($scoreBB{"$idB:$idB"} - $ortoS[$i]);
		}

		# Check if this paralog candidate is already clustered in other clusters
		for $k(1..($i-1)){
		    if (vec($is_paralogB[$k],$hitID,1)){ # Yes, found in cluster $k
			if($debug == 2){
			    print "      $nameB[$hitID] is already in cluster $k, together with:";
			    print " $nameB[$ortoB[$k]] and $nameA[$ortoA[$k]] ";
			    print "($scoreBB{\"$ortoB[$k]:$hitID\"})";
			}
			if (($confPB[$hitID] >= $conf_here) and
			    ($j != 1)){ # The seed ortholog CAN NOT remain there
			    print " and remains there.\n" if ($debug == 2);
			    $paralogs = 0; # No action
			}
			else { # Ortholog of THIS cluster is closer than ortholog of competing cluster $k
			    print " and should be moved here!\n" if ($debug == 2); # Remove from other cluster, add to this cluster
			    @membersBK = split(/ /, $paralogsB[$k]); # This array contains names, not IDs
			    $paralogsB[$k] = "";
			    @tmp = ();
			    for $m(@membersBK){
				push(@tmp,$m) if ($m != $hitID); # Put other members back
			    }
			    $paralogsB[$k] = join(' ',@tmp);
			    undef $is_paralogB[$k]; # Create index that we can check later
			    grep (vec($is_paralogB[$k],$_,1) = 1, @tmp);
			}
			last; # Don't search in other clusters
		    }
		}
		next if (! $paralogs); # Skip that paralog - it is already in cluster $k
		push (@membersB,$hitID);
	    }
	}
	# Calculate confidence values now:
	@tmp = ();
	for $idP (@membersB){ # For each paralog calculate conf value
	    if($scoreBB{"$idB:$idB"} == $ortoS[$i]){
		if ($scoreBB{"$idB:$idP"} == $scoreBB{"$idB:$idB"}){
		    $confPB[$idP] = 1.0;
		}
		else{
		    $confPB[$idP] = 0.0;
		}
	    }
	    else{
		$confPB[$idP] = ($scoreBB{"$idB:$idP"} - $ortoS[$i]) /
		    ($scoreBB{"$idB:$idB"} - $ortoS[$i]);
	    }
	    push (@tmp, $idP) if ($confPB[$idP] >= $conf_cutoff); # If one wishes to use only significant paralogs
	}
	@membersB = @tmp;
	########### Merge if necessary:
	if ($merge[$i] > 0){ # Merge existing cluster with overlapping cluster
	    @tmp = split(/ /,$paralogsB[$merge[$i]]);
	    for $m (@membersB){
		push (@tmp, $m)  unless (vec($is_paralogB[$merge[$i]],$m,1));
	    }
	    $paralogsB[$merge[$i]] = join(' ',@tmp);
	    undef $is_paralogB[$merge[$i]];
	    grep (vec($is_paralogB[$merge[$i]],$_,1) = 1, @tmp); # Refresh index of paralog array
	}
	######### Typical new cluster:
	else{  # Create a new cluster
	    $paralogsB[$i] = join(' ',@membersB);
	    undef $is_paralogB; # Create index that we can check later
	    grep (vec($is_paralogB[$i],$_,1) = 1, @membersB);
	}
    }

#####################################################
    &clean_up(1);
####################################################
# Find group for orphans. If cluster contains only one member, find where it should go:
    for $i (1..$o){
	@membersA = split(/ /, $paralogsA[$i]);
	@membersB = split(/ /, $paralogsB[$i]);
	$na = @membersA;
	$nb = @membersB;
	if (($na == 0) and $nb){
	    print "Warning: empty A cluster $i\n";
	    for $m (@membersB){
		$bestscore = 0;
		$bestgroup = 0;
		$bestmatch = 0;
		for $j (1..$o) {
		    next if ($i == $j); # Really need to check against all 100% members of the group.
		    @membersBJ = split(/ /, $paralogsB[$j]);
		    for $k (@membersBJ){
			next if ($confPB[$k] != 1); # For all 100% in-paralogs
			$score = $scoreBB{"$m:$k"};
			if ($score > $bestscore){
			    $bestscore = $score;
			    $bestgroup = $j;
			    $bestmatch = $k;
			}
		    }
		}
		print "Orphan $nameB[$m] goes to group $bestgroup with $nameB[$bestmatch]\n" ;
		@members = split(/ /, $paralogsB[$bestgroup]);
		push (@members, $m);
		$paralogsB[$bestgroup] = join(' ',@members);
		$paralogsB[$i] = "";
		undef $is_paralogB[$bestgroup];
		undef $is_paralogB[$i];
		grep (vec($is_paralogB[$bestgroup],$_,1) = 1, @members); # Refresh index of paralog array
#		 grep (vec($is_paralogB[$i],$_,1) = 1, ());
	    }
	}
	if ($na and ($nb == 0)){
	    print "Warning: empty B cluster $i\n";
	    for $m (@membersA){
		$bestscore = 0;
		$bestgroup = 0;
		$bestmatch = 0;
		for $j (1..$o) {
		    next if ($i == $j);
		    @membersAJ = split(/ /, $paralogsA[$j]);
		    for $k (@membersAJ){
			next if ($confPA[$k] != 1); # For all 100% in-paralogs
			$score = $scoreAA{"$m:$k"};
			if ($score > $bestscore){
			    $bestscore = $score;
			    $bestgroup = $j;
			    $bestmatch = $k;
			}
		    }
		}
		print "Orphan $nameA[$m] goes to group $bestgroup with $nameA[$bestmatch]\n";
		@members = split(/ /, $paralogsA[$bestgroup]);
		push (@members, $m);
		$paralogsA[$bestgroup] = join(' ',@members);
		$paralogsA[$i] = "";
		undef $is_paralogA[$bestgroup];
		undef $is_paralogA[$i];
		grep (vec($is_paralogA[$bestgroup],$_,1) = 1, @members); # Refresh index of paralog array
#	     grep (vec($is_paralogA[$i],$_,1) = 1, ());
	    }
	}
    }

    &clean_up(1);
###################

# ##############################################################################
# Check for alternative orthologs, sort paralogs by confidence and print results
# ##############################################################################

    for $i(1..$o){
	@membersA = split(/ /, $paralogsA[$i]);
	@membersB = split(/ /, $paralogsB[$i]);
	$message = "";


	$idB = $ortoB[$i];
	$nB = $hitnBA[$idB];
	for $idA(@membersA){
	    next if ($confPA[$idA] != 1.0);
	    $nA = $hitnAB[$idA];
	    $confA[$i] = $ortoS[$i]; # default
	    $bsA[$idA] = 1.0;
	    ##############
	    for $j(1..$nB){
		$idH = $hitBA[$idB][$j];
		################ Some checks for alternative orthologs:
		# 1. Don't consider sequences that are already in this cluster
		next if (vec($is_paralogA[$i],$idH,1));
		next if ($confPA[$idH] > 0); # If $conf_cutoff > 0 idH might be incide circle, but not paralog

		# 2. Check if candidate for alternative ortholog is already clustered in stronger clusters
		$in_other_cluster = 0;
		for $k(1..($i-1)){ # Check if current ortholog is already clustered
		    if (vec($is_paralogA[$k],$idH,1)){
			$in_other_cluster = $k;
			last;
		    }
		}
#		 next if ($in_other_cluster); # This hit is clustered in cluster $k. It cannot be alternative ortholog

		# 3. The best hit of candidate ortholog should be ortoA or at least to belong into this cluster
		@besthits = split (/ /,$besthitAB[$idH]);
		$this_family = 0;
		for $bh (@besthits){
		    $this_family = 1 if ($idB == $bh);
		}
#		 next unless ($this_family); # There was an alternative BA match but it's best match did not belong here
		################# Done with checks - if sequence passed, then it could be an alternative ortholog
		$confA[$i] = $ortoS[$i] - $scoreBA{"$idB:$idH"};

		last;
	    }
	    $message .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.", $nameA[$idA], 100*$bsA[$idA]);
	    $message .= sprintf(" Alternative seed ortholog is %s (%d bits away from this cluster)", $nameA[$idH], $confA[$i]) if ($bsA[$idA] < 0.75);
	    $message .= sprintf("\n");

	}
	########
	$idA = $ortoA[$i];
	$nA = $hitnAB[$idA];
	for $idB(@membersB){
	    next if ($confPB[$idB] != 1.0);
	    $nB = $hitnBA[$idB];
	    $confB[$i] = $ortoS[$i]; # default
	    $bsB[$idB] = 1.0;

	    for $j(1..$nA){ # For all AB hits of given ortholog
		$idH = $hitAB[$idA][$j];
		# ############### Some checks for alternative orthologs:
		# 1. Don't consider sequences that are already in this cluster
		next if (vec($is_paralogB[$i],$idH,1));
		next if ($confPB[$idH] > 0); # If $conf_cutoff > 0 idH might be incide circle, but not paralog

		# 2. Check if candidate for alternative ortholog is already clustered in stronger clusters
		$in_other_cluster = 0;
		for $k(1..($i-1)){
		    if (vec($is_paralogB[$k],$idH,1)){
			$in_other_cluster = $k;
			last; # out from this cycle
		    }
		}
#		 next if ($in_other_cluster); # This hit is clustered in cluster $k. It cannot be alternative ortholog

		# 3. The best hit of candidate ortholog should be ortoA
		@besthits = split (/ /,$besthitBA[$idH]);
		$this_family = 0;
		for $bh (@besthits){
		    $this_family = 1 if ($idA == $bh);
		}
#		 next unless ($this_family); # There was an alternative BA match but it's best match did not belong here
		# ################ Done with checks - if sequence passed, then it could be an alternative ortholog
		$confB[$i] = $ortoS[$i] - $scoreAB{"$idA:$idH"};

		last;
	    }
	    $message .= sprintf("Bootstrap support for %s as seed ortholog is %d%%.", $nameB[$idB], 100*$bsB[$idB]);
	    $message .= sprintf(" Alternative seed ortholog is %s (%d bits away from this cluster)", $nameB[$idH],$confB[$i]) if ($bsB[$idB] < 0.75);
	    $message .= sprintf("\n");


	}
	close FF;


	########### Sort and print members of A ############
	$nA = @membersA;
	$nB = @membersB;
	$nMAX = ($nA > $nB) ? $nA : $nB;
	# Sort membersA inside the cluster by confidence:
	for $m (0..($nA-1)){
	    while($confPA[$membersA[$m]] < $confPA[$membersA[$m+1]]){
		$temp = $membersA[$m];
		$membersA[$m] = $membersA[$m+1];
		$membersA[$m+1] = $temp;
		--$m if ($m > 1);
	    }
	}
	$paralogsA[$i] = join(' ',@membersA); # Put them back together
	# Sort membersB inside the cluster by confidence:
	for $m (0..($nB-1)){
	    while($confPB[$membersB[$m]] < $confPB[$membersB[$m+1]]){
		$temp = $membersB[$m];
		$membersB[$m] = $membersB[$m+1];
		$membersB[$m+1] = $temp;
		--$m if ($m > 1);
	    }
	}
	$paralogsB[$i] = join(' ',@membersB); # Put them back together
	# Print to text file and to HTML file
    }

    if ($table){
	$filename = "$ARGV[2]" . "paranoid_output/" . $ARGV[0] . "." . $ARGV[1] . ".txt";
	open F, ">$filename" or die;
	print F "OrtoID\tScore\tOrtoA\tOrtoB\n";
	for $i(1..$o){
	    print F "$i\t$ortoS[$i]\t";
	    @members = split(/ /, $paralogsA[$i]);
	    for $m (@members){
		$m =~ s/://g;
		printf (F "%s %.3f ", $nameA[$m], $confPA[$m]);
	    }
	    print F "\t";
	    @members = split(/ /, $paralogsB[$i]);
	    for $m (@members){
		$m =~ s/://g;
		printf (F "%s %.3f ", $nameB[$m], $confPB[$m]);
	    }
	    print F "\n";
	}
	close F;
    }

  }

##############################################################
# Functions:
##############################################################
sub clean_up { # Sort members within cluster and clusters by size
############################################################################################### Modification by Isabella 3

    # Sort on index arrays with perl's built in sort instead of using bubble sort.

    $var = shift;
    $totalA = $totalB = 0;
    # First pass: count members within each cluster
    foreach $i (1..$o) {
	@membersA = split(/ /, $paralogsA[$i]);
	$clusnA[$i] = @membersA; # Number of members in this cluster
	$totalA += $clusnA[$i];
	$paralogsA[$i] = join(' ',@membersA);

	@membersB = split(/ /, $paralogsB[$i]);
	$clusnB[$i] = @membersB; # Number of members in this cluster
	$totalB += $clusnB[$i];
	$paralogsB[$i] = join(' ',@membersB);

	$clusn[$i] =  $clusnB[$i] + $clusnA[$i]; # Number of members in given group
    }

    # Create an array used to store the position each element shall have in the final array
    # The elements are initialized with the position numbers
    my @position_index_array = (1..$o);

    # Sort the position list according to cluster size
    my @cluster_sorted_position_list = sort { $clusn[$b] <=> $clusn[$a]} @position_index_array;

    # Create new arrays for the sorted information
    my @new_paralogsA;
    my @new_paralogsB;
    my @new_is_paralogA;
    my @new_is_paralogB;
    my @new_clusn;
    my @new_ortoS;
    my @new_ortoA;
    my @new_ortoB;


   # Add the information to the new arrays in the orer specifeid by the index array
   for (my $index_in_list = 0; $index_in_list < scalar @cluster_sorted_position_list; $index_in_list++) {

	my $old_index = $cluster_sorted_position_list[$index_in_list];

	if (!$clusn[$old_index]) {
		$o = (scalar @new_ortoS) - 1;
		last;
	}

	$new_paralogsA[$index_in_list + 1] = $paralogsA[$old_index];
        $new_paralogsB[$index_in_list + 1] = $paralogsB[$old_index];
    	$new_is_paralogA[$index_in_list + 1] = $is_paralogA[$old_index];
    	$new_is_paralogB[$index_in_list + 1] = $is_paralogB[$old_index];
   	$new_clusn[$index_in_list + 1] = $clusn[$old_index];
	$new_ortoA[$index_in_list + 1] = $ortoA[$old_index];
	$new_ortoB[$index_in_list + 1] = $ortoB[$old_index];
	$new_ortoS[$index_in_list + 1] = $ortoS[$old_index];
   }

    @paralogsA = @new_paralogsA;
    @paralogsB = @new_paralogsB;
    @is_paralogA = @new_is_paralogA;
    @is_paralogB = @new_is_paralogB;
    @clusn = @new_clusn;
    @ortoS = @new_ortoS;
    @ortoA = @new_ortoA;
    @ortoB = @new_ortoB;

    # Create an array used to store the position each element shall have in the final array
    # The elements are initialized with the position numbers
    @position_index_array = (1..$o);

    # Sort the position list according to score
    @score_sorted_position_list = sort { $ortoS[$b] <=> $ortoS[$a] } @position_index_array;

    # Create new arrays for the sorted information
    my @new_paralogsA2 = ();
    my @new_paralogsB2 = ();
    my @new_is_paralogA2 = ();
    my @new_is_paralogB2 = ();
    my @new_clusn2 = ();
    my @new_ortoS2 = ();
    my @new_ortoA2 = ();
    my @new_ortoB2 = ();

   # Add the information to the new arrays in the orer specifeid by the index array
   for (my $index_in_list = 0; $index_in_list < scalar @score_sorted_position_list; $index_in_list++) {

	my $old_index = $score_sorted_position_list[$index_in_list];
	$new_paralogsA2[$index_in_list + 1] = $paralogsA[$old_index];
        $new_paralogsB2[$index_in_list + 1] = $paralogsB[$old_index];
    	$new_is_paralogA2[$index_in_list + 1] = $is_paralogA[$old_index];
    	$new_is_paralogB2[$index_in_list + 1] = $is_paralogB[$old_index];
   	$new_clusn2[$index_in_list + 1] = $clusn[$old_index];
	$new_ortoA2[$index_in_list + 1] = $ortoA[$old_index];
	$new_ortoB2[$index_in_list + 1] = $ortoB[$old_index];
	$new_ortoS2[$index_in_list + 1] = $ortoS[$old_index];
   }

    @paralogsA = @new_paralogsA2;
    @paralogsB = @new_paralogsB2;
    @is_paralogA = @new_is_paralogA2;
    @is_paralogB = @new_is_paralogB2;
    @clusn = @new_clusn2;
    @ortoS = @new_ortoS2;
    @ortoA = @new_ortoA2;
    @ortoB = @new_ortoB2;

#################################################################################### End modification by Isabella 3

}

sub overlap_test {
        my @Fld = @_;

	# Filter out fragmentary hits by:
	# Ignore hit if aggregate matching area covers less than $seq_overlap_cutoff of sequence.
	# Ignore hit if local matching segments cover less than $segment_coverage_cutoff of sequence.
	#
	# $Fld[3] and $Fld[4] are query and subject lengths.
	# $Fld[5] and $Fld[6] are lengths of the aggregate matching region on query and subject. (From start of first matching segment to end of last matching segment).
	# $Fld[7] and $Fld[8] are local matching length on query and subject (Sum of all segments length's on query).

	$retval = 1;
#	if ($Fld[3] >= $Fld[4]) {
		if ($Fld[5] < ($seq_overlap_cutoff * $Fld[3])) {$retval = 0};
		if ($Fld[7] < ($segment_coverage_cutoff * $Fld[3])) {$retval = 0};
#	}

#	if ($Fld[4] >= $Fld[3]) {
		if ($Fld[6] < ($seq_overlap_cutoff * $Fld[4])) {$retval = 0};
		if ($Fld[8] < ($segment_coverage_cutoff * $Fld[4])) {$retval = 0};
#	}

	# print "$Fld[3] $Fld[5] $Fld[7]; $Fld[4] $Fld[6] $Fld[8]; retval=$retval\n";

	return $retval;
}
