#!/usr/bin/perl
open(IN,"hydrogen_bond_summary_5");
open(HOH,">summary_HOH");
open(TYR,">summary_TYR");
open(SER,">summary_SER");
open(THR,">summary_THR");
open(ARG,">summary_ARG");
open(LYS,">summary_LYS");
open(HIS,">summary_HIS");
open(TRP,">summary_TRP");
open(ASN,">summary_ASN");
open(GLN,">summary_GLN");
open(BACK,">summary_backbone");
@in=<IN>;
foreach $line(@in)
{
 if((substr($line,31,3)eq "HOH")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    #$angle2=substr($line,85,8);
    printf HOH ("%4s%8.4f%10.4f\n",$pdb,$dist,$angle1);
   }
 if((substr($line,31,3)eq "TYR")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf TYR ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "SER")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf SER ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "THR")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf THR ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "ARG")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf ARG ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "LYS")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf LYS ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "HIS")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf HIS ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
 if((substr($line,31,3)eq "TRP")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf TRP ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
  if((substr($line,31,3)eq "ASN")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf ASN ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
  if((substr($line,31,3)eq "GLN")&&(substr($line,41,11)ne"backbone NH"))
   {$pdb=substr($line,1,4);
    $dist=substr($line,44,6);
    $angle1=substr($line,67,8);
    $angle2=substr($line,85,8);
    printf GLN ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
  if(substr($line,41,11)eq"backbone NH")
   {$pdb=substr($line,1,4);
    $dist=substr($line,56,6);
    $angle1=substr($line,79,8);
    $angle2=substr($line,97,8);
    printf BACK ("%4s%8.4f%10.4f%10.4f\n",$pdb,$dist,$angle1,$angle2);
   }
} 
