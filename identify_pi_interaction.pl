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
$jj=0;
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
              push(@{$chainID[$i]},$line_3);
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
              push(@{$chainID[$i]},$line_4);
              push(@{$chainID[$i]},$line_c);
              push(@{$chainID[$i]},$line_n);
             } 
           }
        }
     }
  }
print "@{$chainID[$i]}\n";
}
#########identify CN is finished#
@pi=();
$pp=0;
foreach $line(@pdb)
  {if($line=~/^ATOM.{9}(CG|CD1|CD2|CE1|CE2|CZ|NE1|CE3|CZ2|CZ3|CH2|ND1|NE2|NE|NH1|NH2) +(TYR|TRP|PHE|HIS|HIE|HID|ARG)/)
     {$pi[$pp]=$line;
      $pp++;
     }
  }
for($i=0;$i<@chainID;$i++)
{
 @ring=();
 if(@{$chainID[$i]}/3==1)
 { $x_1=substr(${$chainID[$i]}[0],30,8);
   $y_1=substr(${$chainID[$i]}[0],38,8);
   $z_1=substr(${$chainID[$i]}[0],46,8);
   $x_c=substr(${$chainID[$i]}[1],30,8);
   $y_c=substr(${$chainID[$i]}[1],38,8);
   $z_c=substr(${$chainID[$i]}[1],46,8);
   $x_n=substr(${$chainID[$i]}[2],30,8);
   $y_n=substr(${$chainID[$i]}[2],38,8);
   $z_n=substr(${$chainID[$i]}[2],46,8);
   goto lable1;
 }
if((@{$chainID[$i]}/3==2)&&("A" eq "A"))
 { print "The second CN\n"; 
   $x_1=substr(${$chainID[$i]}[3],30,8);
   $y_1=substr(${$chainID[$i]}[3],38,8);
   $z_1=substr(${$chainID[$i]}[3],46,8);
   $x_c=substr(${$chainID[$i]}[4],30,8);
   $y_c=substr(${$chainID[$i]}[4],38,8);
   $z_c=substr(${$chainID[$i]}[4],46,8);
   $x_n=substr(${$chainID[$i]}[5],30,8);
   $y_n=substr(${$chainID[$i]}[5],38,8);
   $z_n=substr(${$chainID[$i]}[5],46,8);
   @ring=(); 
   $vec_1x=$x_1-$x_c;
   $vec_1y=$y_1-$y_c;
   $vec_1z=$z_1-$z_c;##vector from C of CN to attached atom 
   $bond=0;
   $mm=0;
   foreach $line(@lig)
     {if((substr($line,21,1) eq $chainID[$i])&&(substr($line,77,1) ne "H")&&(substr($line,77,1) ne "F")&&(substr($line,76,2) ne "BR")&&(substr($line,76,2) ne "CL")&&(substr($line,76,2) ne " I")&&(substr($line,30,8) != $x_1)&&(substr($line,30,8) != $x_c)&&(substr($line,30,8) != $x_n))
       {$ligand[$mm]=$line;
        $mm++;
       }
     }
   foreach $line(@ligand)
     {
      $x_x=substr($line,30,8);
      $y_x=substr($line,38,8);
      $z_x=substr($line,46,8);
      $distX=sqrt(($x_1-$x_x)**2+($y_1-$y_x)**2+($z_1-$z_x)**2);
      if(($distX<1.8)&&($distX>1.0)) ##bond length standard 
        {
         $atom[$bond]=$line;
         $line=undef;
         $bond++;
        }
     } 
   #print "$bond atom connected to the CCN in chain $chainID[$i]\n";
   if($bond==1)
     {print "Attaced to aphatic CH2 or olefinic CH in chain $chainID[$i]\n";
     }
   elsif($bond==3)
     {print "Attaced to quaternary carbon in chain $chainID[$i]\n";
     }
   elsif($bond==2)
     {$x5=substr($atom[0],30,8);
      $y5=substr($atom[0],38,8);
      $z5=substr($atom[0],46,8);
      $x6=substr($atom[1],30,8);
      $y6=substr($atom[1],38,8);
      $z6=substr($atom[1],46,8);
      $mid_x=($x5+$x6)/2;
      $mid_y=($y5+$y6)/2;
      $mid_z=($z5+$z6)/2;
      $vec_5x=$x5-$x_1;
      $vec_5y=$y5-$y_1;
      $vec_5z=$z5-$z_1;#vector from CN attached atom to atom5
      $vec_6x=$x6-$x_1;
      $vec_6y=$y6-$y_1;
      $vec_6z=$z6-$z_1;#vector from CN attached atom to atom6
      $nomal_x=$vec_5y*$vec_6z-$vec_5z*$vec_6y;
      $nomal_y=$vec_5z*$vec_6x-$vec_5x*$vec_6z;
      $nomal_z=$vec_5x*$vec_6y-$vec_5y*$vec_6x;#nomal vetor of vector5 and vector 6
      $vec_2x=$mid_x-$x_1;
      $vec_2y=$mid_y-$y_1;
      $vec_2z=$mid_z-$z_1;#vector from attached atom to middle point 
      $cos_3=($vec_1x*$nomal_x+$vec_1y*$nomal_y+$vec_1z*$nomal_z)/(sqrt($vec_1x**2+$vec_1y**2+$vec_1z**2)*(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)));
      if(($cos_3<-0.173)&&($cos_3>0.175)) ##the angle between $vec_1 and $nomal is less than 80,more than 100
        {print "CN is Attached to tertiary carbon CH in chain $chainID[$i]\n";
        }
      elsif(($cos_3>-0.173)&&($cos_3<0.175))##the angle between $vec_1 and $nomal is less than 100 more than 80
        {$aa=0;
         $bb=0;
         foreach $line_2(@ligand)
           {
            $x_y=substr($line_2,30,8);
            $y_y=substr($line_2,38,8);
            $z_y=substr($line_2,46,8);
            $dist5=sqrt(($x5-$x_y)**2+($y5-$y_y)**2+($z5-$z_y)**2);
            $dist6=sqrt(($x6-$x_y)**2+($y6-$y_y)**2+($z6-$z_y)**2);
            if(($dist5<1.8)&&($dist5>1.20))##double bond standard
              {$atom5[$aa]=$line_2;
               $aa++;
               $line_2=undef;
              }
            if(($dist6<1.8)&&($dist6>1.20))##double bond standard
              {$atom6[$bb]=$line_2;
               $bb++;
               $line_2=undef;
              }
           }
         #print "Attached number is $aa and $bb\n";
         if($aa==0||$bb==0)
           {print "CN is Attached to a olefinic carbon C in chain $chainID[$i]\n";
           }
         if($aa==3||$bb==3)
           {print "It is not a ring\n";
           }
         if(($aa==1)&&($bb==1))
          {
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);#atom51 is attached to atom5
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5; #vector from atom51 to atom5
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);#atom61 is attached to atom6
           $vec_4x=$x61-$x6;
           $vec_4y=$y61-$y6;
           $vec_4z=$z61-$z6;#vector from atom61 to atom6
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $cos_4=($vec_3x*$vec_4x+$vec_3y*$vec_4y+$vec_3z*$vec_4z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));##the angle 
           $cos_5=($nomal_x*$vec_3x+$nomal_y*$vec_3y+$nomal_z*$vec_3z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)));#the angle between plane normal vector and vector4
           $cos_6=($nomal_x*$vec_4x+$nomal_y*$vec_4y+$nomal_z*$vec_4z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));#the same as above
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175)) #the angle more than 80less than 100,the distance is double bond
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 111\n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];
             }
           if(($dist_7>2.25)&&($dist_7<2.62)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))#the $dist_7 is the distance skip a atom in ring
             {foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_8=sqrt(($x51-$x_z)**2+($y51-$y_z)**2+($z51-$z_z)**2);
                 $dist_9=sqrt(($x61-$x_z)**2+($y61-$y_z)**2+($z61-$z_z)**2);
                 if(($dist_8<1.80)&&($dist_8>1.20)&&($dist_9<1.80)&&($dist_9>1.20))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 111\n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;
                   }
                }
             }
           }
         if(($aa==1)&&($bb==2))
           {@tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $vec_4x=$x61-$x6;
           $vec_4y=$y61-$y6;
           $vec_4z=$z61-$z6;
           $x62=substr($atom6[1],30,8);
           $y62=substr($atom6[1],38,8);
           $z62=substr($atom6[1],46,8);
           $vec_6x=$x62-$x6;
           $vec_6y=$y62-$y6;
           $vec_6z=$z62-$z6;
           $cos_5=($nomal_x*$vec_3x+$nomal_y*$vec_3y+$nomal_z*$vec_3z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)));
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x51-$x62)**2+($y51-$y62)**2+($z51-$z62)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_4x+$vec_3y*$vec_4y+$vec_3z*$vec_4z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));
              $cos_6=($nomal_x*$vec_4x+$nomal_y*$vec_4y+$nomal_z*$vec_4z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));
               foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);#@dist is array which 
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 222 \n";

                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;


                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_3x*$vec_6x+$vec_3y*$vec_6y+$vec_3z*$vec_6z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)));
              $cos_6=($nomal_x*$vec_6x+$nomal_y*$vec_6y+$nomal_z*$vec_6z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)));
               foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 222 \n";

                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175)) 
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 222 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];
             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 222 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[1];
             }

          }
        if(($aa==2)&&($bb==2))
          {
           @tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $x52=substr($atom5[1],30,8);
           $y52=substr($atom5[1],38,8);
           $z52=substr($atom5[1],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $vec_4x=$x52-$x5;
           $vec_4y=$y52-$y5;
           $vec_4z=$z52-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $x62=substr($atom6[1],30,8);
           $y62=substr($atom6[1],38,8);
           $z62=substr($atom6[1],46,8);
           $vec_5x=$x61-$x6;
           $vec_5y=$y61-$y6;
           $vec_5z=$z61-$z6;
           $vec_6x=$x62-$x6;
           $vec_6y=$y62-$y6;
           $vec_6z=$z62-$z6;
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x51-$x62)**2+($y51-$y62)**2+($z51-$z62)**2);
           $dist_9=sqrt(($x52-$x61)**2+($y52-$y61)**2+($z52-$z61)**2);
           $dist_10=sqrt(($x52-$x62)**2+($y52-$y62)**2+($z52-$z62)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8,$dist_9,$dist_10;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_5x+$vec_3y*$vec_5y+$vec_3z*$vec_5z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));      
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_12=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;
                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_3x*$vec_6x+$vec_3y*$vec_6y+$vec_3z*$vec_6z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_6x*$nomal_x+$vec_6y*$nomal_y+$vec_6z*$nomal_z)/(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_12=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             } 
           if($tmp[0]==$dist_9)
             {$cos_4=($vec_4x*$vec_5x+$vec_4y*$vec_5y+$vec_4z*$vec_5z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_12=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   } 
                }
             }
           if($tmp[0]==$dist_10)
             {$cos_4=($vec_4x*$vec_6x+$vec_4y*$vec_6y+$vec_4z*$vec_6z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_6x*$nomal_x+$vec_6y*$nomal_y+$vec_6z*$nomal_z)/(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_12=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];

             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[1];

             }
           if(($dist_9<1.80)&&($dist_9>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[0];

             }
           if(($dist_10<1.80)&&($dist_10>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[1];

             }

          }
        if(($aa==2)&&($bb==1))
          {@tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $x52=substr($atom5[1],30,8);
           $y52=substr($atom5[1],38,8);
           $z52=substr($atom5[1],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $vec_4x=$x52-$x5;
           $vec_4y=$y52-$y5;
           $vec_4z=$z52-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $vec_5x=$x61-$x6;
           $vec_5y=$y61-$y6;
           $vec_5z=$z61-$z6;
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x52-$x61)**2+($y52-$y61)**2+($z52-$z61)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_5x+$vec_3y*$vec_5y+$vec_3z*$vec_5z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 444 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_4x*$vec_5x+$vec_4y*$vec_5y+$vec_4z*$vec_5z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 444 \n";
                    $ring[0]=${$chainID[$i]}[3];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 444 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];

             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 444 \n";
              $ring[0]=${$chainID[$i]}[3];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[0];

             }

          }
        }
     }
  
#######identify pi interaction######
if(@ring)
  {#print "@ring\n";
   $r_x=0;
   $r_y=0;
   $r_z=0;
   foreach $line(@ring)
     {$r_x=$r_x+substr($line,30,8);
      $r_y=$r_y+substr($line,38,8);
      $r_z=$r_z+substr($line,46,8);
     }
   $r_x=$r_x/@ring;
   $r_y=$r_y/@ring;
   $r_z=$r_z/@ring;##center of the ring
   $v_r1x=substr($ring[2],30,8)-substr($ring[0],30,8);
   $v_r1y=substr($ring[2],38,8)-substr($ring[0],38,8);
   $v_r1z=substr($ring[2],46,8)-substr($ring[0],46,8);
   $v_r2x=substr($ring[4],30,8)-substr($ring[2],30,8);
   $v_r2y=substr($ring[4],38,8)-substr($ring[2],38,8);
   $v_r2z=substr($ring[4],46,8)-substr($ring[2],46,8);
   $nomal_rx=$v_r1y*$v_r2z-$v_r1z*$v_r2y;
   $nomal_ry=$v_r1z*$v_r2x-$v_r1x*$v_r2z;
   $nomal_rz=$v_r1x*$v_r2y-$v_r1y*$v_r2x;##nomal of the ring
   $res_now=substr($pdb[0],17,9);
   foreach $line(@pi)
     {$res=substr($line,17,9);
      if($res eq $res_now)
        {push(@array,$line);
        }
      else
        {
         $p_x=0;
         $p_y=0;
         $p_z=0;
         if(@array>4)
         {foreach $line(@array)
           {$p_x=$p_x+substr($line,30,8);
            $p_y=$p_y+substr($line,38,8);
            $p_z=$p_z+substr($line,46,8);
           }
          
          $p_x=$p_x/@array;
          $p_y=$p_y/@array;
          $p_z=$p_z/@array;
          $dist_pi=sqrt(($r_x-$p_x)**2+($r_y-$p_y)**2+($r_z-$p_z)**2);
          $v_p1x=substr($array[2],30,8)-substr($array[0],30,8);
          $v_p1y=substr($array[2],38,8)-substr($array[0],38,8);
          $v_p1z=substr($array[2],46,8)-substr($array[0],46,8);
          $v_p2x=substr($array[4],30,8)-substr($array[2],30,8);
          $v_p2y=substr($array[4],38,8)-substr($array[2],38,8);
          $v_p2z=substr($array[4],46,8)-substr($array[2],46,8);
          $nomal_px=$v_p1y*$v_p2z-$v_p1z*$v_p2y;
          $nomal_py=$v_p1z*$v_p2x-$v_p1x*$v_p2z;
          $nomal_pz=$v_p1x*$v_p2y-$v_p1y*$v_p2x;##the nomal of ring residues
          $cos_pi=($nomal_rx*$nomal_px+$nomal_ry*$nomal_py+$nomal_rz*$nomal_pz)/(sqrt($nomal_rx**2+$nomal_ry**2+$nomal_rz**2)*sqrt($nomal_px**2+$nomal_py**2+$nomal_pz**2));
          $angle=acos($cos_pi)*180/pi;
          if(($dist_pi<5.0)&&((($cos_pi>0.93)&&($cos_pi<1))||(($cos_pi>-1)&&($cos_pi<-0.93))))
            {print "The ring parallel pi interaction with $res_now in chain $chainID[$i] with distance $dist_pi and angle $angle \n";  ###sprintf("%.2f", $number)
            }
          if(($dist_pi<5.5)&&(($cos_pi>-0.34)&&($cos_pi<0.34)))
            {print "The ring edge to face pi interaction with $res_now in chain $chainID[$i] with distance $dist_pi and angle $angle \n";
            }
          }
          $res_now=$res;
          @array=();
          push(@array,$line);
        }
     }
  }
 }
if(("A" ne "C")&&(@{$chainID[$i]}/3==2))
 {print "The first CN\n";
  $x_1=substr(${$chainID[$i]}[0],30,8);
  $y_1=substr(${$chainID[$i]}[0],38,8);
  $z_1=substr(${$chainID[$i]}[0],46,8);
  $x_c=substr(${$chainID[$i]}[1],30,8);
  $y_c=substr(${$chainID[$i]}[1],38,8);
  $z_c=substr(${$chainID[$i]}[1],46,8);
  $x_n=substr(${$chainID[$i]}[2],30,8);
  $y_n=substr(${$chainID[$i]}[2],38,8);
  $z_n=substr(${$chainID[$i]}[2],46,8);
  goto lable1;
 }

lable1:
   @ring=();
   $vec_1x=$x_1-$x_c;
   $vec_1y=$y_1-$y_c;
   $vec_1z=$z_1-$z_c;##vector from C of CN to attached atom 
   $bond=0;
   $mm=0;
   foreach $line(@lig)
     {if((substr($line,21,1) eq $chainID[$i])&&(substr($line,77,1) ne "H")&&(substr($line,77,1) ne "F")&&(substr($line,76,2) ne "BR")&&(substr($line,76,2) ne "CL")&&(substr($line,76,2) ne " I")&&(substr($line,30,8) != $x_1)&&(substr($line,30,8) != $x_c)&&(substr($line,30,8) != $x_n))
       {$ligand[$mm]=$line;
        $mm++;
       }
     }
   foreach $line(@ligand)
     {
      $x_x=substr($line,30,8);
      $y_x=substr($line,38,8);
      $z_x=substr($line,46,8);
      $distX=sqrt(($x_1-$x_x)**2+($y_1-$y_x)**2+($z_1-$z_x)**2);
      if(($distX<1.8)&&($distX>1.0)) ##bond length standard 
        {
         $atom[$bond]=$line;
         $line=undef;
         $bond++;
        }
     } 
   #print "$bond atom connected to the CCN in chain $chainID[$i]\n";
   if($bond==1)
     {print "Attaced to aphatic CH2 or olefinic CH in chain $chainID[$i]\n";
     }
   elsif($bond==3)
     {print "Attaced to quaternary carbon in chain $chainID[$i]\n";
     }
   elsif($bond==2)
     {$x5=substr($atom[0],30,8);
      $y5=substr($atom[0],38,8);
      $z5=substr($atom[0],46,8);
      $x6=substr($atom[1],30,8);
      $y6=substr($atom[1],38,8);
      $z6=substr($atom[1],46,8);
      $mid_x=($x5+$x6)/2;
      $mid_y=($y5+$y6)/2;
      $mid_z=($z5+$z6)/2;
      $vec_5x=$x5-$x_1;
      $vec_5y=$y5-$y_1;
      $vec_5z=$z5-$z_1;#vector from CN attached atom to atom5
      $vec_6x=$x6-$x_1;
      $vec_6y=$y6-$y_1;
      $vec_6z=$z6-$z_1;#vector from CN attached atom to atom6
      $nomal_x=$vec_5y*$vec_6z-$vec_5z*$vec_6y;
      $nomal_y=$vec_5z*$vec_6x-$vec_5x*$vec_6z;
      $nomal_z=$vec_5x*$vec_6y-$vec_5y*$vec_6x;#nomal vetor of vector5 and vector 6
      $vec_2x=$mid_x-$x_1;
      $vec_2y=$mid_y-$y_1;
      $vec_2z=$mid_z-$z_1;#vector from attached atom to middle point 
      $cos_3=($vec_1x*$nomal_x+$vec_1y*$nomal_y+$vec_1z*$nomal_z)/(sqrt($vec_1x**2+$vec_1y**2+$vec_1z**2)*(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)));
      if(($cos_3<-0.173)&&($cos_3>0.175)) ##the angle between $vec_1 and $nomal is less than 80,more than 100
        {print "CN is Attached to tertiary carbon CH in chain $chainID[$i]\n";
        }
      elsif(($cos_3>-0.173)&&($cos_3<0.175))##the angle between $vec_1 and $nomal is less than 100 more than 80
        {$aa=0;
         $bb=0;
         foreach $line_2(@ligand)
           {
            $x_y=substr($line_2,30,8);
            $y_y=substr($line_2,38,8);
            $z_y=substr($line_2,46,8);
            $dist5=sqrt(($x5-$x_y)**2+($y5-$y_y)**2+($z5-$z_y)**2);
            $dist6=sqrt(($x6-$x_y)**2+($y6-$y_y)**2+($z6-$z_y)**2);
            if(($dist5<1.8)&&($dist5>1.20))##double bond standard
              {$atom5[$aa]=$line_2;
               $aa++;
               $line_2=undef;
              }
            if(($dist6<1.8)&&($dist6>1.20))##double bond standard
              {$atom6[$bb]=$line_2;
               $bb++;
               $line_2=undef;
              }
           }
         #print "Attached number is $aa and $bb\n";
         if($aa==0||$bb==0)
           {print "CN is Attached to a olefinic carbon C in chain $chainID[$i]\n";
           }
         if($aa==3||$bb==3)
           {print "It is not a ring\n";
           }
         if(($aa==1)&&($bb==1))
          {
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);#atom51 is attached to atom5
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5; #vector from atom51 to atom5
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);#atom61 is attached to atom6
           $vec_4x=$x61-$x6;
           $vec_4y=$y61-$y6;
           $vec_4z=$z61-$z6;#vector from atom61 to atom6
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $cos_4=($vec_3x*$vec_4x+$vec_3y*$vec_4y+$vec_3z*$vec_4z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));##the angle 
           $cos_5=($nomal_x*$vec_3x+$nomal_y*$vec_3y+$nomal_z*$vec_3z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)));#the angle between plane normal vector and vector4
           $cos_6=($nomal_x*$vec_4x+$nomal_y*$vec_4y+$nomal_z*$vec_4z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));#the same as above
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175)) #the angle more than 80less than 100,the distance is double bond
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 111\n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];
             }
           if(($dist_7>2.25)&&($dist_7<2.62)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))#the $dist_7 is the distance skip a atom in ring
             {foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_8=sqrt(($x51-$x_z)**2+($y51-$y_z)**2+($z51-$z_z)**2);
                 $dist_9=sqrt(($x61-$x_z)**2+($y61-$y_z)**2+($z61-$z_z)**2);
                 if(($dist_8<1.80)&&($dist_8>1.20)&&($dist_9<1.80)&&($dist_9>1.20))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 111\n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;
                   }
                }
             }
           }
         if(($aa==1)&&($bb==2))
           {@tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $vec_4x=$x61-$x6;
           $vec_4y=$y61-$y6;
           $vec_4z=$z61-$z6;
           $x62=substr($atom6[1],30,8);
           $y62=substr($atom6[1],38,8);
           $z62=substr($atom6[1],46,8);
           $vec_6x=$x62-$x6;
           $vec_6y=$y62-$y6;
           $vec_6z=$z62-$z6;
           $cos_5=($nomal_x*$vec_3x+$nomal_y*$vec_3y+$nomal_z*$vec_3z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)));
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x51-$x62)**2+($y51-$y62)**2+($z51-$z62)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_4x+$vec_3y*$vec_4y+$vec_3z*$vec_4z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));
              $cos_6=($nomal_x*$vec_4x+$nomal_y*$vec_4y+$nomal_z*$vec_4z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)));
               foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);#@dist is array which 
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 222 \n";

                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;


                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_3x*$vec_6x+$vec_3y*$vec_6y+$vec_3z*$vec_6z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)));
              $cos_6=($nomal_x*$vec_6x+$nomal_y*$vec_6y+$nomal_z*$vec_6z)/(sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2)*(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)));
               foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 222 \n";

                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175)) 
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 222 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];
             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 222 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[1];
             }

          }
        if(($aa==2)&&($bb==2))
          {
           @tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $x52=substr($atom5[1],30,8);
           $y52=substr($atom5[1],38,8);
           $z52=substr($atom5[1],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $vec_4x=$x52-$x5;
           $vec_4y=$y52-$y5;
           $vec_4z=$z52-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $x62=substr($atom6[1],30,8);
           $y62=substr($atom6[1],38,8);
           $z62=substr($atom6[1],46,8);
           $vec_5x=$x61-$x6;
           $vec_5y=$y61-$y6;
           $vec_5z=$z61-$z6;
           $vec_6x=$x62-$x6;
           $vec_6y=$y62-$y6;
           $vec_6z=$z62-$z6;
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x51-$x62)**2+($y51-$y62)**2+($z51-$z62)**2);
           $dist_9=sqrt(($x52-$x61)**2+($y52-$y61)**2+($z52-$z61)**2);
           $dist_10=sqrt(($x52-$x62)**2+($y52-$y62)**2+($z52-$z62)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8,$dist_9,$dist_10;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_5x+$vec_3y*$vec_5y+$vec_3z*$vec_5z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));      
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_12=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;
                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_3x*$vec_6x+$vec_3y*$vec_6y+$vec_3z*$vec_6z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_6x*$nomal_x+$vec_6y*$nomal_y+$vec_6z*$nomal_z)/(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_12=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             } 
           if($tmp[0]==$dist_9)
             {$cos_4=($vec_4x*$vec_5x+$vec_4y*$vec_5y+$vec_4z*$vec_5z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_12=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   } 
                }
             }
           if($tmp[0]==$dist_10)
             {$cos_4=($vec_4x*$vec_6x+$vec_4y*$vec_6y+$vec_4z*$vec_6z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_6x*$nomal_x+$vec_6y*$nomal_y+$vec_6z*$nomal_z)/(sqrt($vec_6x**2+$vec_6y**2+$vec_6z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_11=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_12=sqrt(($x_z-$x62)**2+($y_z-$y62)**2+($z_z-$z62)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_11>1.2)&&($dist_11<1.80)&&($dist_12>1.2)&&($dist_12<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 333 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[1];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];

             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[1];

             }
           if(($dist_9<1.80)&&($dist_9>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[0];

             }
           if(($dist_10<1.80)&&($dist_10>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 333 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[1];

             }

          }
        if(($aa==2)&&($bb==1))
          {@tmp=();
           $x51=substr($atom5[0],30,8);
           $y51=substr($atom5[0],38,8);
           $z51=substr($atom5[0],46,8);
           $x52=substr($atom5[1],30,8);
           $y52=substr($atom5[1],38,8);
           $z52=substr($atom5[1],46,8);
           $vec_3x=$x51-$x5;
           $vec_3y=$y51-$y5;
           $vec_3z=$z51-$z5;
           $vec_4x=$x52-$x5;
           $vec_4y=$y52-$y5;
           $vec_4z=$z52-$z5;
           $x61=substr($atom6[0],30,8);
           $y61=substr($atom6[0],38,8);
           $z61=substr($atom6[0],46,8);
           $vec_5x=$x61-$x6;
           $vec_5y=$y61-$y6;
           $vec_5z=$z61-$z6;
           $dist_7=sqrt(($x51-$x61)**2+($y51-$y61)**2+($z51-$z61)**2);
           $dist_8=sqrt(($x52-$x61)**2+($y52-$y61)**2+($z52-$z61)**2);
           @tmp=sort{$a<=>$b}$dist_7,$dist_8;
           if($tmp[0]==$dist_7)
             {$cos_4=($vec_3x*$vec_5x+$vec_3y*$vec_5y+$vec_3z*$vec_5z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_3x*$nomal_x+$vec_3y*$nomal_y+$vec_3z*$nomal_z)/(sqrt($vec_3x**2+$vec_3y**2+$vec_3z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x51)**2+($y_z-$y51)**2+($z_z-$z51)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 444 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[0];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   }
                }
             }
           if($tmp[0]==$dist_8)
             {$cos_4=($vec_4x*$vec_5x+$vec_4y*$vec_5y+$vec_4z*$vec_5z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2));
              $cos_5=($vec_4x*$nomal_x+$vec_4y*$nomal_y+$vec_4z*$nomal_z)/(sqrt($vec_4x**2+$vec_4y**2+$vec_4z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              $cos_6=($vec_5x*$nomal_x+$vec_5y*$nomal_y+$vec_5z*$nomal_z)/(sqrt($vec_5x**2+$vec_5y**2+$vec_5z**2)*sqrt($nomal_x**2+$nomal_y**2+$nomal_z**2));
              foreach $line_3(@ligand)
                {$x_z=substr($line_3,30,8);
                 $y_z=substr($line_3,38,8);
                 $z_z=substr($line_3,46,8);
                 $dist_9=sqrt(($x_z-$x52)**2+($y_z-$y52)**2+($z_z-$z52)**2);
                 $dist_10=sqrt(($x_z-$x61)**2+($y_z-$y61)**2+($z_z-$z61)**2);
                 if(($tmp[0]>2.25)&&($tmp[0]<2.62)&&($dist_9>1.2)&&($dist_9<1.80)&&($dist_10>1.2)&&($dist_10<1.80)&&($cos_4>0.94)&&($cos_4<1)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
                   {print "CN is Attached to a six plane ring in chain $chainID[$i] 444 \n";
                    $ring[0]=${$chainID[$i]}[0];
                    $ring[1]=$atom[0];
                    $ring[2]=$atom[1];
                    $ring[3]=$atom5[1];
                    $ring[4]=$atom6[0];
                    $ring[5]=$line_3;

                   }
                }
             }
           if(($dist_7<1.80)&&($dist_7>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 444 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[0];
              $ring[4]=$atom6[0];

             }
           if(($dist_8<1.80)&&($dist_8>1.20)&&($cos_5>-0.173)&&($cos_5<0.175)&&($cos_6>-0.173)&&($cos_6<0.175))
             {print "CN is Attached to a five plane ring in chain $chainID[$i] 444 \n";
              $ring[0]=${$chainID[$i]}[0];
              $ring[1]=$atom[0];
              $ring[2]=$atom[1];
              $ring[3]=$atom5[1];
              $ring[4]=$atom6[0];

             }

          }
        }
     }
  
#######identify pi interaction######
if(@ring)
  {#print "@ring\n";
   $r_x=0;
   $r_y=0;
   $r_z=0;
   foreach $line(@ring)
     {$r_x=$r_x+substr($line,30,8);
      $r_y=$r_y+substr($line,38,8);
      $r_z=$r_z+substr($line,46,8);
     }
   $r_x=$r_x/@ring;
   $r_y=$r_y/@ring;
   $r_z=$r_z/@ring;##center of the ring
   $v_r1x=substr($ring[2],30,8)-substr($ring[0],30,8);
   $v_r1y=substr($ring[2],38,8)-substr($ring[0],38,8);
   $v_r1z=substr($ring[2],46,8)-substr($ring[0],46,8);
   $v_r2x=substr($ring[4],30,8)-substr($ring[2],30,8);
   $v_r2y=substr($ring[4],38,8)-substr($ring[2],38,8);
   $v_r2z=substr($ring[4],46,8)-substr($ring[2],46,8);
   $nomal_rx=$v_r1y*$v_r2z-$v_r1z*$v_r2y;
   $nomal_ry=$v_r1z*$v_r2x-$v_r1x*$v_r2z;
   $nomal_rz=$v_r1x*$v_r2y-$v_r1y*$v_r2x;##nomal of the ring
   $res_now=substr($pdb[0],17,9);
   foreach $line(@pi)
     {$res=substr($line,17,9);
      if($res eq $res_now)
        {push(@array,$line);
        }
      else
        {
         $p_x=0;
         $p_y=0;
         $p_z=0;
         if(@array>4)
         {foreach $line(@array)
           {$p_x=$p_x+substr($line,30,8);
            $p_y=$p_y+substr($line,38,8);
            $p_z=$p_z+substr($line,46,8);
           }
          
          $p_x=$p_x/@array;
          $p_y=$p_y/@array;
          $p_z=$p_z/@array;
          $dist_pi=sqrt(($r_x-$p_x)**2+($r_y-$p_y)**2+($r_z-$p_z)**2);
          $v_p1x=substr($array[2],30,8)-substr($array[0],30,8);
          $v_p1y=substr($array[2],38,8)-substr($array[0],38,8);
          $v_p1z=substr($array[2],46,8)-substr($array[0],46,8);
          $v_p2x=substr($array[4],30,8)-substr($array[2],30,8);
          $v_p2y=substr($array[4],38,8)-substr($array[2],38,8);
          $v_p2z=substr($array[4],46,8)-substr($array[2],46,8);
          $nomal_px=$v_p1y*$v_p2z-$v_p1z*$v_p2y;
          $nomal_py=$v_p1z*$v_p2x-$v_p1x*$v_p2z;
          $nomal_pz=$v_p1x*$v_p2y-$v_p1y*$v_p2x;##the nomal of ring residues
          $cos_pi=($nomal_rx*$nomal_px+$nomal_ry*$nomal_py+$nomal_rz*$nomal_pz)/(sqrt($nomal_rx**2+$nomal_ry**2+$nomal_rz**2)*sqrt($nomal_px**2+$nomal_py**2+$nomal_pz**2));
          $angle=acos($cos_pi)*180/pi;
          if(($dist_pi<5.0)&&((($cos_pi>0.93)&&($cos_pi<1))||(($cos_pi>-1)&&($cos_pi<-0.93))))
            {print "The ring parallel pi interaction with $res_now in chain $chainID[$i] with distance $dist_pi and angle $angle \n";  ###sprintf("%.2f", $number)
            }
          if(($dist_pi<5.5)&&(($cos_pi>-0.34)&&($cos_pi<0.34)))
            {print "The ring edge to face pi interaction with $res_now in chain $chainID[$i] with distance $dist_pi and angle $angle \n";
            }
          }
          $res_now=$res;
          @array=();
          push(@array,$line);
        }
     }
  }
}
