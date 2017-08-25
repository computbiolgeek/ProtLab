#!/usr/bin/perl

## global variables
$pdb_dir = "/dors/meilerlab/apps/PISCES/opm/2014/all/";
$bcl_executable = "bcl.exe PDBConvert ";
$out_dir = "/dors/meilerlab/apps/PISCES/opm/2014/";
$global_pdb_dir = "/dors/csb/pdb/distr/";

mkdir $out_dir;
mkdir $pdb_dir;

## check the parameters and print usage if not
if($#ARGV != 0)
{
  die "usage: <pdb_list>\n"
      ."<pdb_list> : list of 5 letter codes\n";
}

my @files;
if( -f $ARGV[0])
{
	## get the list of pdb ids from the given file
	open( IN, $ARGV[0]) or die "cannot open pdb list at $ARGV[0]\n";
	@files=<IN>;
	close IN;
}
else
{
	@files=@ARGV;
}

$membrane_string  = "  <MEMBRANE>\n";
$membrane_string .= "    <NORMAL X=\"0.00\" Y=\"0.00\" Z=\"####\"/>\n";
$membrane_string .= "    <TMATRIX>\n";
$membrane_string .= "      <ROWX X=\"1.00\" Y=\"0.00\" Z=\"0.00\" T=\"0.00\"/>\n";
$membrane_string .= "      <ROWY X=\"0.00\" Y=\"1.00\" Z=\"0.00\" T=\"0.00\"/>\n";
$membrane_string .= "      <ROWZ X=\"0.00\" Y=\"0.00\" Z=\"1.00\" T=\"0.00\"/>\n";
$membrane_string .= "    </TMATRIX>\n";
$membrane_string .= "  </MEMBRANE>\n";

## for every file
foreach $file_id( @files)
{
  ## chomp and lower case
  chomp $file_id;
  print $file_id."\n";

  ## ensure correct length is provided
  if( length( $file_id) != 5){ die "The provided pdb code should be 5 letters not $file_id\n";}
  
  ## parse pdb and chain id and mid_id
  $pdb_id = substr( $file_id, 0, 4);
  $pdb_id = lc( $pdb_id);
  $chain_id = substr( $file_id, 4, 1);
  $mid_id = lc( substr( $pdb_id, 1, 2));
  
  $bcl_chain_id = $chain_id;
  if( $chain_id =~ /_/)
  {
    $bcl_chain_id = "' '";
  }

  ## now generate the fasta file
  ## first get the pdb filename
  $pdb_filename = $pdb_dir.$pdb_id.".pdb";
  $local_out_dir = $out_dir.$mid_id."/";
  if( ! -f $pdb_filename)
  {
    print "PDB with code: ".$pdb_id." is not present in ".$pdb_dir.". Attempting to download from OPM...";
    system( "curl -m 60 -s http://opm.phar.umich.edu/pdb/".$pdb_id.".pdb -o ".$pdb_filename);
    open my $pdbfile, '<', "$pdb_filename"; 
    my $firstLine = <$pdbfile>; 
    close $pdbfile;
    if( $firstLine =~ /^REMARK/)
    {
      print "Success\n";
    }
    elsif( $firstLine =~ /<html>/)
    {
      print "Failed. Goto http://opm.phar.umich.edu/server.php and save the pdb\n";
      unlink $pdb_filename;
      next;
    }
  } 
 
  $xml_file = "$local_out_dir".$pdb_id.".xml"; 
  if( ! -f $xml_file)
  {
    # write out fake xml file 
    open my $pdbfile, '<', "$pdb_filename"; 
    my $firstLine = <$pdbfile>; 
    my @firstLineTok = split(/:/,$firstLine);
    chomp $firstLineTok[-1];
    my $membrane_thickness = $firstLineTok[-1]+0.0;
    close $pdbfile;
    $membrane_string_copy = $membrane_string;
    $membrane_string_copy =~ s/####/$membrane_thickness/;
    open my $pdbfile, '>', "$xml_file"; 
    print $pdbfile $membrane_string_copy."\n";
    close $pdbfile;
  }

  ## assemble command and execute it and store the error
  mkdir $local_out_dir;
  $local_pdb_filename=$local_out_dir.$pdb_id.".pdb";
  $log_file=$local_out_dir.$pdb_id.".bcl.log";
  if( ! -f $local_pdb_filename)
  {
  	open my $pdbfile, '<', "$pdb_filename";
	    $global_pdbfile=$global_pdb_dir.$mid_id."/pdb".$pdb_id.".ent.gz";
	    if( ! -f $global_pdbfile)
	    {
	      system( "cat ".$pdb_filename." | grep '^REMARK' > ".$local_pdb_filename);
              print "PDB Distribution filename: ".$global_pdbfile." does not exist! This file needs to be present to obtain the SEQRES information that OPM omits or gets wrong. Better hope its right in this case.\n";
              system( "cat ".$pdb_filename." | egrep -v '^(REMARK|ANISOU)' | egrep -v '^HETATM...........(HOH|DUM)' >> ".$local_pdb_filename);
	    }
	    else
	    {
              system( "cat ".$pdb_filename." | grep '^REMARK' | grep -v '^REMARK 350' > ".$local_pdb_filename);
               #$get_chains_opm = "grep '^ATOM ' $pdb_filename | awk '{if( ".'$3 == "N"){print substr($0,22,1)}}\' | sort | uniq | tr -d \'\\n\'';
	       #$get_chains_pdb = "zgrep '^ATOM ' $global_pdbfile | awk '{if( ".'$3 == "N"){print substr($0,22,1)}}\' | sort | uniq | tr -d \'\\n\'';
	       #$chains_opm = qx( $get_chains_opm);
	       #$chains_pdb = qx( $get_chains_pdb);
	       #if($chains_opm ne $chains_pdb)
	       #{
               #  # print $pdb_id." has chains ".$chains_opm." in the opm version but has chains ".$chains_pdb." in the pdb version.\n";
	       #  system( "echo ".$pdb_id." chains in opm: ".$chains_opm." chains in pdb: ".$chains_pdb);
	       #  unlink($local_pdb_filename);
               #  exit;
	       #}
	      system( "bcl.exe PDBConvert ".$global_pdbfile." -bcl_pdb -biomolecule -output_prefix " .$local_out_dir.$pdb_id.".generate_chains. -logger File ".$log_file);
	      system( "grep '^SEQRES' " .$local_out_dir.$pdb_id.".generate_chains.bcl.pdb >> ".$local_pdb_filename);
              # remove the generated file, to avoid naming collisions and to avoid polluting the file system
              unlink $local_out_dir.$pdb_id.".generate_chains.bcl.pdb";
	      system( "cat ".$pdb_filename." | egrep -v '^(REMARK|ANISOU|SEQRES)' | egrep -v '^HETATM...........(HOH|DUM)' >> ".$local_pdb_filename);
            }
  }
  if( ! -f $local_out_dir.$pdb_id.$chain_id.".pdb")
  {
    $command = $bcl_executable.$local_pdb_filename." -aaclass AAComplete -chains ".$bcl_chain_id." -output_prefix ".$local_out_dir.$pdb_id." -bcl_pdb Split -logger File ".$log_file;
    $error = system( $command);

    ## issue warning about error
    if( $error)
    {
       print "pdb reading failed: $file_id with error $error using command: $command\n";
    }
  }
  $fasta_file = "$local_out_dir".$pdb_id.$chain_id.".fasta"; 
  if( ! -f $fasta_file)
  {
    ## assemble command and execute it and store the error
    $command = $bcl_executable.$local_pdb_filename." -fasta -chains ".$bcl_chain_id." -bcl_pdb Split -output_prefix ".$local_out_dir.$pdb_id." -convert_to_natural_aa_type -aaclass AAComplete -logger File ".$log_file;
    $error = system( $command);

    ## issue warning about error
    if( $error)
    {
      print "pdb reading failed: $file_id with error $error using command: $command\n";
    }
  }
  
}  

# print "Now run the following command (preferrably on silicon or somewhere with a lot of spare processors";
# print "find ".$out_dir." -name '????.pdb' | xargs -n1 -P24 /dors/meilerlab/apps/scripts/sequence_analysis/SecondaryStructureAnalysis.py -no-analysis -multimer";


