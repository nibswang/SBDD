#!/usr/bin/perl
open(PDB,"$ARGV[0]");
open(LIST,"list_lig_new");
open(OUT,">>list_identify_CN");
@pdb=<PDB>;
@list=<LIST>;
$j=0;
foreach $line(@pdb)
  {$i=0;
   foreach $line2(@list)
     {$id=substr($line2,0,3);
      if ($line=~/^HETATM...........$id/)
        {$i++;
        }
     }
   if($i==1)
     {#print "$line";
      $lig[$j]=$line;
      $j++;
     }
  }
$n=0;
$nn=0;
$chainID[$n]=substr($lig[0],21,1);
foreach $line(@lig)
  {
   $chain=substr($line,21,1);
   foreach $line_5(@chainID) 
     {if ($chain ne $line_5)
        {$n++;
         push(@chainID,$chain);
        }
     }
  }
$l=0;
$m=0;
foreach $line(@lig)
  {
   if((substr($line,13,1) eq "C")&&(substr($line,21,1)eq $chainID))
     {
      $carbon[$l]=$line;
      $l++;
     }
   elsif((substr($line,13,1) eq "N")&&(substr($line,21,1)eq $chainID))
     {
      $nitrogen[$m]=$line;
      $m++;
     }
  }
foreach $line_n(@nitrogen)
  {
   $atom_n=substr($line_n,11,5);
   $xn=substr($line_n,30,8);
   $yn=substr($line_n,38,8);
   $zn=substr($line_n,46,8);
   foreach $line_c(@carbon)
     {$atom_c=substr($line_c,11,5);
      $xc=substr($line_c,30,8);
      $yc=substr($line_c,38,8);
      $zc=substr($line_c,46,8);
      $distC_N=sqrt(($xn-$xc)**2+($yn-$yc)**2+($zn-$zc)**2);
      $vet_CN_x=$xc-$xn;
      $vet_CN_y=$yc-$yn;
      $vet_CN_z=$zc-$zn;##vetor from N to C
      if (($distC_N>1.0)&&($distC_N<1.35))
        {
         $number=0;
         foreach $line_3(@carbon)
           {$cos_0=0;
            $x3=substr($line_3,30,8);
            $y3=substr($line_3,38,8);
            $z3=substr($line_3,46,8);
            $distC_C=sqrt(($xc-$x3)**2+($yc-$y3)**2+($zc-$z3)**2);
           ####angle judge####
            $vet_CC_x=$x3-$xc;
            $vet_CC_y=$y3-$yc;
            $vet_CC_z=$z3-$zc;##vector from C of CN to C
            $xxx=(sqrt($vet_CN_x**2+$vet_CN_y**2+$vet_CN_z**2))*(sqrt($vet_CC_x**2+$vet_CC_y**2+$vet_CC_z**2));
            if ($xxx!=0)
              {
                $cos_0=($vet_CN_x*$vet_CC_x+$vet_CN_y*$vet_CC_y+$vet_CN_z*$vet_CC_z)/$xxx;
              }
            if (($distC_C>1.2)&&($distC_C<1.60)&&($cos_0<1)&&($cos_0>0.8660))  ##bond length and bond angle 
             {
              $number++;
             }
           }
         foreach $line_4(@nitrogen)
           {$cos_1=0;
            $x4=substr($line_4,30,8);
            $y4=substr($line_4,38,8);
            $z4=substr($line_4,46,8);
            if($x4!=$xn)
              {$distN_C=sqrt(($xc-$x4)**2+($yc-$y4)**2+($zc-$z4)**2);
               $vet_NC_x=$x4-$xc;
               $vet_NC_y=$y4-$yc;
               $vet_NC_z=$z4-$zc;
               $cos_1=($vet_CN_x*$vet_NC_x+$vet_CN_y*$vet_NC_y+$vet_CN_z*$vet_NC_z)/((sqrt($vet_CN_x**2+$vet_CN_y**2+$vet_CN_z**2))*(sqrt($vet_NC_x**2+$vet_NC_y**2+$vet_NC_z**2)));
              }
           if(($distN_C>1.2)&&($distN_C<1.4)&&($cos_1<1)&&($cos_1>0.866))
             {
              $number++;
             } 
            
           }
         if ($number==1)
           { print "$ARGV[0]:Nitrile\n";
             $cn[0]=$line_n;
             $cn[1]=$line_c;
             #print "@cn";
             print OUT ("$ARGV[0]\n");
           }
        }
     }
  }




























































































