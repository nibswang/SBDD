#!/usr/bin/perl
use Math::Trig; 
open(PDB,"$ARGV[0]");
open(LIST,"list_lig_new");
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
#####chain###
$nn=0;
$chainID[0]=substr($lig[0],21,1);
foreach $line(@lig)
  {$n=0;
   $chain=substr($line,21,1);
   foreach $line_5(@chainID)
     {if ($chain eq $line_5)
        {
         $n++;
        }
     }
   if($n==0)
     {$nn++;
      $chainID[$nn]=$chain;
     }
  }
################
for($i=0;$i<@chainID;$i++)
{
$l=0;
$m=0;
@carbon=();
@nitrogen=();
foreach $line(@lig)
  {
   if((substr($line,13,1) eq "C")&&(substr($line,21,1)eq $chainID[$i]))
     {
      $carbon[$l]=$line;
      $l++;
     }
   elsif((substr($line,13,1) eq "N")&&(substr($line,21,1)eq $chainID[$i]))
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
      if (($distC_N>0.9)&&($distC_N<1.45))
        {
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
            if (($distC_C>1.1)&&($distC_C<1.70)&&($cos_0<1)&&($cos_0>0.8660))  ##bond length and bond angle 
             {
              print "$ARGV[0]:$chainID[$i]:Nitrile\n";
             # push(@{$chainID[$i]},$line_3);
              push(@{$chainID[$i]},$line_c);
              push(@{$chainID[$i]},$line_n);
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
           if(($distN_C>1.1)&&($distN_C<1.6)&&($cos_1<1)&&($cos_1>0.866))
             {
              print "$ARGV[0]:$chainID[$i]:Nitrile\n";
             # push(@{$chainID[$i]},$line_4);
              push(@{$chainID[$i]},$line_c);
	      push(@{$chainID[$i]},$line_n);
	     } 
           }
        }
     }
  }
}
#########identify CN is finished#
$jj=0;
#open(OUT,">out.pdb");
foreach $line(@pdb)
  {if($line=~/^ATOM/)
        {$rec[$jj]=$line;
         $jj++;
        }
        
  }
############################
for($i=0;$i<@chainID;$i++)
  {for($k=0;$k< @{$chainID[$i]}/2;$k++)
     {$x_c[$k]=substr(${$chainID[$i]}[$k*2],30,8);
      $y_c[$k]=substr(${$chainID[$i]}[$k*2],38,8);
      $z_c[$k]=substr(${$chainID[$i]}[$k*2],46,8);
      $x_n[$k]=substr(${$chainID[$i]}[$k*2+1],30,8);
      $y_n[$k]=substr(${$chainID[$i]}[$k*2+1],38,8);
      $z_n[$k]=substr(${$chainID[$i]}[$k*2+1],46,8);
      foreach $line(@rec)
        {$x=substr($line,30,8);
         $y=substr($line,38,8);
         $z=substr($line,46,8);
         $dist_hy2=sqrt(($x-$x_n[$k])**2+($y-$y_n[$k])**2+($z-$z_n[$k])**2);
         if($dist_hy2<4.0)
           {
            $hash{$dist_hy2}=$line;
           }
        }
       @key=keys %hash;
       @tmp=sort{$a<=>$b}@key;
       $atom_type=substr($hash{$tmp[0]},13,1);
      if($atom_type eq "C")
           {print "The $k Nitrile form hydrophobic interaction in chain $chainID[$i]\n";
            print "$tmp[0]\n";
            print "$hash{$tmp[0]}";
           }
        
     }
  }
