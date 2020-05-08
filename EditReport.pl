#!/usr/bin/perl

## Author: Dhawal Jain (Park lab, HMS)
## The script reads the hyperedited reads and reports edited position with respect to the reference genome 

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use constant { true => 1, false => 0 };
BEGIN { our $start_run = time(); }

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $bam = "";
my $hd ="";
my $type ="";
my $wd = ""; 
my $out ="";
my $help = 0;
Getopt::Long::GetOptions(
  'bam=s'            => \$bam,
  'hd:s'             => \$hd,
  'type:s'           => \$type,
  'out:s'            => \$out,
  'wd=s'             => \$wd,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if ($j) {
   print "\nUsage: perl EditReport.pl -bam [FILE_PATH] -hd [FILE_PATH] -out [OUTPREFIX] -type [AG or TC] -wd [PATH]\n\n";
   print "This script parses output of the mapped hyperedited sites\n\n";
   print "Options (required*):\n";
   print "   -bam             Input bam file\n";
   print "   -hd              Header file for output bam\n";
   print "   -out             outprefix\n";   
   print "   -type            Mutation type (AG or TC)\n";   
   print "\n";
   print "Options (optional):\n";
   print "   -wd              Directory to write the output file into\n";
   print "\n";
   print "   -help|-h         Display usage information.\n";
   print "\n";
   print "Default outputs:\n";
   print "    Writes a bam and txt file\n\n\n";
   exit 0;
  }
}
if($help or $bam ne "-" or $hd eq "" or $out eq "" or $type eq ""){
  print " One or more inputs are not recognized**\n";
  help(1);
}
$wd =~ s/\/$//;

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print  "[EditReport unmapped reads] START:\t $run_time seconds\n";
print  " Command: perl EditReport.pl -bam $bam -hd $hd -out $out -wd $wd\n";

my $baseout = $out;
$out = $wd."/".$out if($wd ne "");

my $outl = $out.".report.txt.gz";
open(OUTL,"| gzip -c - > $outl") or die $!;

my $outa = $out.".anal.txt.gz";
open(OUTA,"| gzip -c - > $outa") or die $!;

my $out2 = $out.".vis.bam";
open O2, "| samtools view -bS - > $out2" or die "can't create $out2";
if($hd ne ""){
  open(F, "$hd") or die $!;
  my @header = <F>;
  print O2 @header;
  close F;
}

my $FIN="STDIN";
if($bam eq "-"){
   print "Input is pipein \n";
   $FIN="STDIN";
}elsif($bam ne "-" and $bam =~ m/.bam$/){
  open($FIN,"samtools view -@ 4 -h $bam | ") or die "no bam: $bam";
  print " Input bam file: $bam\n";
}

my %cntr;
while(<$FIN>){
  next if(/^(\#)/);
  next if(/^(\@)/); 
  chomp;
  s/\r//;  
  my @temp = split(/\t/);
  next if($temp[1]>=256);
  $cntr{"tot"}++;
  next if($temp[1] ==4);
  
  $cntr{"mapq0"}++;
  $cntr{"mapq10"}++ if($temp[4]>=10);
  $cntr{"mapq60"}++ if($temp[4]>=60);
  next if($temp[4]<60);
  
  if($cntr{"tot"}%1000000 ==0){
    print " processed $cntr{tot} lines.\n";
  }
  
  my ($md) = $_ =~ /\tMD:Z:(\d*)/;
  my ($nm) = "0";
  if($_=~ /\tNM:i:(\d*)/){ $nm=$1;}
  
  my ($aa) = $temp[5]=~ /^(\d+)S\S+/;
  $aa=0 if(!defined($aa) || $aa eq"");  ## left hand side
  
  my $mm =mismatch($md,$nm);
  my $end = get_read_mapping_span($temp[3],$temp[5]);
  editAnalysis($temp[0],$temp[2],$temp[3],$temp[5],$temp[4],$temp[9],$mm,$aa, $end, $type); 
  #$id,$chr,$start,$cigar,$mapq,$seq, mm, aclip, end, type
  
  my $id= $temp[0];
  $id =~ s/^.*\|//;
  $id = "OP:Z:".$id;
  $temp[0]=~ s/\|.*//;
  push @temp, $id;
  print O2 join("\t",@temp),"\n";

}

print " Sample\tTotalMates\tMapped(>=0)\tMapped(>=10)\tMapped(>=60)\n";
print " $baseout\t$cntr{tot}\t$cntr{mapq0}\t$cntr{mapq10}\t$cntr{mapq60}\n";
$watch_run = time();
$run_time = $watch_run - $start_run;
print "[EditReport] END:\t $run_time seconds\n";

close($FIN);
close (O2);
close(OUTL);
exit 0;


###-------------------------------------------------------------------------------------------------------------------------------
# subroutines
###-------------------------------------------------------------------------------------------------------------------------------
sub editAnalysis {
  my ($id,$chr,$start,$cigar,$mapq,$seq, $mm, $aa,$end, $type) =@_;
 
  $id=~ s/^.*\|//;
  my @rde = split(";",$id );
  
  if($type eq "TC"){
    foreach my $i (0..$#rde){
      $rde[$i] = (length($seq) - $rde[$i] -1);
    }
    @rde = reverse @rde;
  }
  
  ## remove alterations that are actually mismatches  
  if($mm ne ""){
    my @mm = split(";",$mm);
    if($aa>0){
      foreach my $i(0..$#mm){ $mm[$i] += ($aa-1);  }
    }
    my %mm = map { $_ => 1 } @mm;  
    foreach my $i (0..$#rde){
      if($mm{$rde[$i]}){
         delete $rde[$i];
      }
    }
    @rde = grep defined, @rde; 
  }
    
  my @len = split (/\D+/,$cigar); # storing the length per operation
  my @ops = split (/\d+/,$cigar); # storing the operation
  shift @ops; # remove the empty first element
  my $ops_string="";
  foreach my $i (1 .. scalar@ops){
    if($ops[$i-1] eq "M" or $ops[$i-1] eq "D"){
      $ops_string .=$ops[$i-1] x $len[$i-1];
    }else{
      $ops_string .="-" x $len[$i-1];
    }
  }
  my @ops_string = split("",$ops_string);
  my @genomicpos;
  my $k=0;
  foreach my $i (0..$#ops_string){
     if($ops_string[$i] eq "M" ){
       push @genomicpos, ($start+$k);
       $k++;
     }elsif($ops_string[$i] eq "D" ){
       push @genomicpos, "-";      
       $k++;
     }
     else{
       push @genomicpos, "-";      
     } 
  }

  if(scalar@rde>=4){
    my $iid = $chr.$start.$cigar;
    my @ts;
    foreach my $i (@rde){
      if($genomicpos[$i] ne "-"){
        my $pos = $genomicpos[$i];
        push @ts, $pos;
        #print "$chr\t$pos\t$pos\t$iid\t$mapq\t+\n";
        print OUTL "$chr\t$pos\t$pos\t$iid\t$mapq\t+\n";
      }
    }
    my $l = scalar@ts;
    my $ts = join(",",@ts);
    #print "$chr\t$start\t$end\t$iid\t$mapq\t$type\t$l\n";
    print OUTA "$chr\t$start\t$end\t$iid\t$mapq\t$type\t$l\n";
  }
}
  
sub mismatch {
  my ($md,$nm) =@_;
  return("") if($nm==0);
  my @ops = split ("([\\^]*[ACGT]+)[0]*",$md); # storing the length per operation
  my $sq="";
  foreach my $i(@ops){
    if($i ne "" and is_int($i) ){
      $sq .= "-" x $i;
    }elsif($i =~ m/\^/){
      $sq .= "-" x (length($i)-1);
    }elsif(length($i)==1){
      $sq .= "M";
    }
  }
  return(patternMatch($sq,"M"));
}

sub patternMatch {
  my ($seq,$motif) = @_;
  my @locs;
  while ($seq =~ /$motif/g){
    push @locs, $-[0];
  }
  my $res = "";
  $res = join(";",@locs) if(scalar @locs >0);
  return($res)
}

sub is_int { 
    my $str = $_[0]; 
    $str =~ s/^\s+|\s+$//g;          
    if ($str =~ /^(\-|\+)?\d+?$/) {
        return(1);
    }
    else{
        return(0);
    }
}  

sub get_read_mapping_span {
  my ($start,$cigar) =@_;
  if($cigar eq "*"){
    return($start);
  }

  my @len1 = $cigar =~ m/(\d+)D/g;
  my @len2 = $cigar =~ m/(\d+)M/g;
  
  my $length=0;
  foreach (@len1) { $length += $_; }
  foreach (@len2) { $length += $_; }
  return(($length+$start));
}
