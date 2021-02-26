#!/usr/bin/perl
open(PDB,"$ARGV[0]");
open(LIST,"list_error_2");
open(ST,">>list_simple_PDB");
@pdb=<PDB>;
@list=<LIST>;
$i=0;
foreach $line(@pdb)
  {# $res=substr($line,17,3);
   # print "The residue is $res\n";
    foreach $line2(@list)
     {
      $id=substr($line2,0,3);
     # print "The list is $id\n";
      if($line=~/HETATM...........$id/)
        {print "$ARGV[0]\n";
         print "$line";
         printf ST ("$ARGV[0]\n");
         $i++;
        }
     }
   if ($i!=0)
     {last;
     }
  }

