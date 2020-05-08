#!/usr/bin/perl
# Author:
# Dhawal Jain, Dept of Biomedical Informatics, Harvard

# This program merges clusters on opposite strands and identifies Duplicated Target Site (TSD) from the supported reads

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use constant { true => 1, false => 0 };
use utf8;

my $seq ="GTGGGGGGGTGCCGGCGGGGCTGGGTGGTGGCGGCCGTGTGGCTGGTGGGCTTTTTGTCGTCTCGGGGGTCCCGGGGGGCCTTCTGGTGGGGGGCTGGGGCCCGGGCCGGGGGCGCGGGGGGGGCTGCTGCTGTGGGGCCTTGGGGGGCC";
my $cigar = "141M9S";#"114M36S";
my $start = 3216830;
my $id = "A00873:21:HCWC2DRXX:2:2256:18195:10520_2|5;12;13;14;28;32;33;35;37;40;41;45;46;53;58;61;72;74;75;82;90;104;109;113;117;119;122;125;126;130;131;133;136;141;142;143;145;146";
my $chr = " chr1";
my $mapq = 60;
my $md = "14C5T9T20T0C4G2C34T11G2";
my $nm = "9";

  my ($aa) = $cigar=~ /^(\d+)S\S+/;
  $aa=0 if(!defined($aa) || $aa eq"");  ## left hand side
  
  my $mm =mismatch($md,$nm);
  my $end = get_read_mapping_span($start,$cigar);
  editAnalysis($id,$chr,$start,$cigar,$mapq,$seq,$mm,$aa, $end, "TC"); 
  


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
        print "$chr\t$pos\t$pos\t$iid\t$mapq\t+\n";
        #print OUTL "$chr\t$pos\t$pos\t$iid\t$mapq\t+\n";
      }
    }
    my $l = scalar@ts;
    my $ts = join(",",@ts);
    print "$chr\t$start\t$end\t$iid\t$mapq\t$type\t$l\n";
    #print OUTA "$chr\t$start\t$end\t$iid\t$mapq\t$type\t$l\n";
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
