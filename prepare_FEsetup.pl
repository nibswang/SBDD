#!/usr/bin/perl
## prepare files for FEP calculations, including topology files, position restraint files


open(TOP,"morph.top");
@top=<TOP>;
foreach $line(@top)
{@tmp=split(/\s+/,$line);
 if((@tmp==9)&&($tmp[3] eq "1")&&($tmp[1]== $tmp[6])&&(($tmp[-1]==14.0100)||($tmp[-1]==1.0080)||($tmp[-1]==12.0100)||($tmp[-1]==32.0600)||($tmp[-1]==12.0100)||($tmp[-1]==16.0000)))
   {push @atom,$line;
   }
}

$seq_res=0;
for($i=0;$i<@atom;$i++)
{@tmp=split(/\s+/,$atom[$i]); 
 if($tmp[5] eq "N")
   {
    for($j=$i;$j<$i+26;$j++)
      {@tmp1=split(/\s+/,$atom[$j]);
       if($tmp1[5]eq "O")
         {
          @tmp2=split(/\s+/,$atom[$j+1]);
          if($tmp2[5] eq "OXT")
            {
             $seq_res++;
             for($k=$i;$k<=$j+1;$k++)
                {
                 #print "$atom[$k]";
                 @tmp3=split(/\s+/,$atom[$k]);
                 $tmp3[3]=$seq_res;
                 #printf "%7d%4s%7d%6s%5s%8d%12.6f%12.6f\n",$tmp3[1],$tmp3[2],$tmp3[3],$tmp3[4],$tmp3[5],$tmp3[6],$tmp3[7],$tmp3[8];
                 $new_line=sprintf("%7d%4s%7d%6s%5s%8d%12.6f%12.6f\n",$tmp3[1],$tmp3[2],$tmp3[3],$tmp3[4],$tmp3[5],$tmp3[6],$tmp3[7],$tmp3[8]);
                 #print $new_line;
                 push @new_atom,$new_line;
                }
             #print "\n";
             last;
            }
          else
            {
             $seq_res++;
             for($k=$i;$k<=$j;$k++)
               {
                #print "$atom[$k]";
                @tmp3=split(/\s+/,$atom[$k]);
                $tmp3[3]=$seq_res;
                #printf "%7d%4s%7d%6s%5s%8d%12.6f%12.6f\n",$tmp3[1],$tmp3[2],$tmp3[3],$tmp3[4],$tmp3[5],$tmp3[6],$tmp3[7],$tmp3[8];
                $new_line=sprintf("%7d%4s%7d%6s%5s%8d%12.6f%12.6f\n",$tmp3[1],$tmp3[2],$tmp3[3],$tmp3[4],$tmp3[5],$tmp3[6],$tmp3[7],$tmp3[8]);
                #print $new_line;
                push @new_atom,$new_line;
               }
             #print "\n";
             last;
            }
         }
       
      }
   }
}
#print "@new_atom";

open(NEW,">new_morph.top");
for($i=0;$i<@top;$i++)
{@tmp=split(/\s+/,$top[$i]);
 if(($top[$i]=~/nr  type resno resnm atom    cgnr      charge        mass      typeB    chargeB      massB/)&&($top[$i+1]!~/Cl-/)&&($top[$i+1]!~/Na+/))
   {print "add new atoms here!\n";
    print NEW ($top[$i]);
    print NEW (@new_atom);
   }
 elsif((@tmp==9)&&($tmp[3] eq "1")&&($tmp[1]== $tmp[6])&&(($tmp[-1]==14.0100)||($tmp[-1]==1.0080)||($tmp[-1]==12.0100)||($tmp[-1]==32.0600)||($tmp[-1]==12.0100)||($tmp[-1]==16.0000)))
   {
    $aa=0;
   }
 elsif(($top[$i]=~/ifdef POSRES/)&&($top[$i+1]=~/posres_MOL1.itp/)&&($top[$i+2]=~/endif/))
   {print "MOL1.itp\n";
    print NEW ($top[$i+1]);
    $i=$i+3;
   }
 else
   {print NEW ($top[$i]);
   }
}
close NEW;


####new topology for complex is done##


open(NOLIG,"ligand_removed.pdb");
@nolig=<NOLIG>;
close NOLIG;

open(POS,">posres_MOL1.itp");
print POS ("[position_restraints]\n");

foreach $line(@nolig)
{
 @tmp=split(/\s+/,$line);
 if((@tmp==9)&&($tmp[3]ne "WAT")&&($tmp[3]ne "WAT")&&($tmp[3]ne "Na+")&&($tmp[3]ne "Cl-")&&($tmp[-1] ne "H"))
   {#print "@tmp\n";
    printf POS "%7d%4d%6d%6d%6d\n",$tmp[1],1,1000,1000,1000;
   }
}
close NOLIG;
####position restraint file for protein is done###

open(ITP,"pert.itp");
@itp=<ITP>;
close ITP;

open(NEWTOP,">new_pert.itp");
open(POSL,">posres_LIG.itp");
print POSL ("[position_restraints]\n");


foreach $line(@itp)
{
 print NEWTOP ($line);
 @tmp=split(/\s+/,$line);
 if((@tmp==12)&&($tmp[4]eq "LIG")&&($tmp[8] ne "1.0080"))
   {
    printf POSL "%7d%4d%6d%6d%6d\n",$tmp[1],1,1000,1000,1000;
   }
}
print NEWTOP ("\#include \"posres_LIG\.itp\"");
close NEWTOP;
close POSL;

###position restraint file for ligand is done, and the position restrint was added to ligand topology file###

`mv morph.top old_morph.top`;
`mv new_morph.top morph.top`;

`mv pert.itp old_pert.itp`;
`mv new_pert.itp pert.itp`;

`cp morph.top prod_morph.top`;
`sed -i '/posres_MOL1.itp/d' prod_morph.top`;
`sed -i 's/pert.itp/prod_pert.itp/g' prod_morph.top`;

`cp pert.itp prod_pert.itp`;
`sed -i '/posres_LIG.itp/d' prod_pert.itp `;

`gmx editconf -f morph.gro -o restraint.gro`;
