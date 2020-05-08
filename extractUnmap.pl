#!/usr/bin/perl

## Author: Dhawal Jain (Park lab, HMS)
## The script extracts unmapped mates from the paired end bam file
 ## edits the read ID and the sequence (A>G)
## output is a fastq file ready for remapping  

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
my $wd = ""; 
my $help = 0;
Getopt::Long::GetOptions(
  'bam=s'            => \$bam,
  'wd=s'             => \$wd,
  'help'             => \$help,
  'h'                => \$help,
) or die "Incorrect input! Use -h for usage.\n";
sub help{
  my $j = shift;
  if ($j) {
   print "\nUsage: perl extractUnmap.pl -bam [FILE_PATH] -wd [PATH]\n\n";
   print "This script extracts unmapped reads from bam file, edit the reads and write them as fastq file\n\n";
   print "Options (required*):\n";
   print "   -bam              Input bam file\n";
   print "                       (Input can be piped, else reading bam file will require samtools in the path variable)";
   print "\n";
   print "Options (optional):\n";
   print "   -wd                Directory to write the output file into\n";
   print "\n";
   print "   -help|-h          Display usage information.\n";
   print "\n";
   print "Default outputs:\n";
   print "    Writes fastq file, gzipped format \n\n\n";
   exit 0;
  }
}
help($help);
help(!$bam);
$wd =~ s/\/$//;

###-------------------------------------------------------------------------------------------------------------------------------
# inputs etc
###-------------------------------------------------------------------------------------------------------------------------------
my $outf = $bam;
$outf =~ s/.*\///g;
$outf =~ s/.bam/.fq.gz/;
$outf = $wd."/".$outf if($wd ne "");
open(OUT,"| gzip -c - > $outf") or die $!;

my $outl = $bam;
$outl =~ s/.bam/.logs/;
$outl =~ s/.*\///g;
$outl = $wd."/".$outl if($wd ne "");
open(OUTL,"> $outl") or die $!;

my $watch_run = time();
my $run_time = $watch_run - our $start_run;
print OUTL "[extract unmapped reads] START:\t $run_time seconds\n";
print OUTL " Command: perl extractUnmap.pl -bam $bam -wd $wd\n";

if($bam=~ /.bam/){
  open IN,"samtools view -@ 6 -f 4 $bam|" or next "Can't open bamfile $bam";
}else{
  open IN,"cat $bam|" or next "Can't open file $bam";
}

my %cntr;
my $ct=0;
while(<IN>){
  next if(/^(\#)/);
  next if(/^(\@)/); 
  chomp;
  s/\r//;  
  my @temp = split(/\t/);
  next if($temp[1]>=256);

  if($temp[1] & 0x40){
     $temp[0] .="_1";
     $cntr{"read"}++;
  }elsif($temp[1] & 0x80){
     $temp[0] .="_2";
     $cntr{"mate"}++;
  }else{
    next;
  }
  $ct++;
  if($ct%1000000 ==0){
    print OUTL " processed $ct lines.\n";
  }
  my $poly = polymerCheck($temp[9],25);
  next if($poly eq "polyA" or $poly eq "polyT");
  
  #my $pm = patternMatch($temp[9],"A");
  #$temp[9] =~ s/A/G/g;
  
  my $pm = patternMatch($temp[9],"T");
  $temp[9] =~ s/T/C/g;
  $temp[0] .="|".$pm;

  print OUT "@"."$temp[0]\n$temp[9]\n+\n$temp[10]\n"; 
  #print "@"."$temp[0]\n$temp[9]\n+\n$temp[10]\n"; 

}

print OUTL " Number of reads in pairs: $cntr{read}\n";
print OUTL " Number of mates in pairs: $cntr{mate}\n";
my $x = ($cntr{read} + $cntr{mate});
print OUTL " Total reads processed: $x\n";
$watch_run = time();
$run_time = $watch_run - $start_run;
print OUTL "[extract unmapped reads] END:\t $run_time seconds\n";


close(IN);
close(OUTL);
close(OUT);
exit 0;


###-------------------------------------------------------------------------------------------------------------------------------
# subroutines
###-------------------------------------------------------------------------------------------------------------------------------
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

sub polymerCheck{
  my ($seq,$cut) =@_;
  my $start = 0;
  my $end = length($seq);
  my %fq;
  $fq{A} = $seq =~ tr/A|a//;
  $fq{T} = $seq =~ tr/T|t//;
  $fq{G} = $seq =~ tr/G|g//;
  $fq{C} = $seq =~ tr/C|c//;
  my @keys = keys %fq;
  @keys = sort{ $fq{$b} <=> $fq{$a} or $b cmp $a} @keys if(scalar @keys >1);  
  my$rat = round( ($fq{$keys[0]}*100)/length($seq),2);
  
  my $out ="*";
  if($rat>=$cut){
     #$out = "poly".$keys[0]."(".$rat."%)";
     $out = "poly".$keys[0];
  }
  return($out);
}

sub round {
  my ($n, $places) = @_;
  my $abs = abs $n;
  my $val = substr($abs + ('0.' . '0' x $places . '5'),
                   0,
                   length(int($abs)) +
                     (($places > 0) ? $places + 1 : 0)
                  );
  ($n < 0) ? "-" . $val : $val;
}
