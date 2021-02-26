#!/usr/bin/perl
use Math::Trig;
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
      if (($distC_N>0.9)&&($distC_N<1.45))##nitrile bond length
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
            if (($distC_C>1.1)&&($distC_C<1.7)&&($cos_0<1)&&($cos_0>0.8660))  ##bond length and bond angle 
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
           if(($distN_C>1.1)&&($distN_C<1.6)&&($cos_1<1)&&($cos_1>0.866)) ##bond length and angle
             {
              $number++;
             } 
            
           }
         if ($number==1)
           { print "$ARGV[0]:$chainID[$i]:Nitrile\n";
             #$cn[$j]=$line_n;
	     #$cn[$j+1]=$line_c;
	     #$j=$j+2;
	     #@{$chainID[$i]}=($line_n,$line_c);
	     push(@{$chainID[$i]},$line_n);
             push(@{$chainID[$i]},$line_c);
             #print "@{$chainID[$i]}";#every array contain N and C of nitrile
             #print "${$chainID[$i]}[0]";
             print OUT ("$ARGV[0]\n");
           }
        }
     }
  }
}
##########
$jj=0;
foreach $line(@pdb)
{
 if ($line=~/^ATOM/)
  {
   $rec[$jj]=$line;
   $jj++;
  }
 elsif(($line=~/^HETATM...........HOH/)||($line=~/^HETATM...........SOL/))
  {$rec[$jj]=$line;
   $jj++;
  }
}
  
for($i=0;$i<@chainID;$i++)
{
 for($k=0;$k<@{$chainID[$i]}/2;$k++)
 {
 $x_n[$k]=substr(${$chainID[$i]}[$k*2],30,8);
 $y_n[$k]=substr(${$chainID[$i]}[$k*2],38,8);
 $z_n[$k]=substr(${$chainID[$i]}[$k*2],46,8);
 $x_c[$k]=substr(${$chainID[$i]}[$k*2+1],30,8);
 $y_c[$K]=substr(${$chainID[$i]}[$k*2+1],38,8);
 $z_c[$k]=substr(${$chainID[$i]}[$k*2+1],46,8);
 $vector1x=$x_c[$k]-$x_n[$k];
 $vector1y=$y_c[$K]-$y_n[$k];
 $vector1z=$z_c[$k]-$z_n[$k]; 
 #foreach $line(@rec)
 for($z=0;$z<@rec;$z++)
   {$res=substr($rec[$z],17,9);
    $x_r=substr($rec[$z],30,8);
    $y_r=substr($rec[$z],38,8);
    $z_r=substr($rec[$z],46,8);
    #print "It is $res \n";
    $dist_N=sprintf("%.4f",sqrt(($x_r-$x_n[$k])**2+($y_r-$y_n[$k])**2+($z_r-$z_n[$k])**2));
    $dist_C=sprintf("%.4f",sqrt(($x_r-$x_c[$k])**2+($y_r-$y_c[$k])**2+($z_r-$z_c[$k])**2));
    $angle1=0;
    $angle2=0;
    if(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "OH ")) #TYR
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2)); 
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
          {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CZ "))
             {$connect_x=substr($rec[$y],30,8);
              $connect_y=substr($rec[$y],38,8);
              $connect_z=substr($rec[$y],46,8);
              $vector3x=$connect_x-$x_r;
              $vector3y=$connect_y-$y_r;
              $vector3z=$connect_z-$z_r;
             }
          } 
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "OG "))#SER
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3) eq "CB "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "OG1")) #THR
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CB "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "N  ")) ##backbone NH
      {$res_num=substr($rec[$z],22,4);
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];

       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {$res_num2=substr($rec[$y],22,4);
          $rec_x=substr($rec[$y],30,8);
          $rec_y=substr($rec[$y],38,8);
          $rec_z=substr($rec[$y],46,8);
          $dist_amide=sqrt(($rec_x-$x_r)**2+($rec_y-$y_r)**2+($rec_z-$z_r)**2);
          if(($dist_amide>1.28)&&($dist_amide<1.38)&&(substr($rec[$y],13,3)eq "C  "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s backbone NH:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res backbone NH:$dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&((substr($rec[$z],13,3) eq "NH1")||(substr($rec[$z],13,3) eq "NH2")||(substr($$rec[$z],13,3) eq "NE "))) #ARG
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CZ "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "NZ "))#LYS
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CE "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "ND2")) # ASN
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CG "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res : $dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "NE2")&&(substr($rec[$z],17,3) eq "GLN")) # GLN
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CD "))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2)); ##439
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res : $dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }

    elsif(($dist_N<3.4)&&((substr($rec[$z],13,3) eq "NE2")||(substr($rec[$z],13,3) eq "ND1"))&&((substr($rec[$z],17,3) eq "HIS")||(substr($rec[$z],17,3) eq "HIE")||(substr($rec[$z],17,3) eq "HID")))#HIS
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CE1"))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res : $dist_N A, with angle1 $angle1,angle2 $angle2\n";

      }
    elsif(($dist_N<3.4)&&(substr($rec[$z],13,3) eq "NE1")) #TRP
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       for($y=$z-75;$y<$z+75;$y++)
         {if((substr($rec[$y],17,9)eq $res)&&(substr($rec[$y],13,3)eq "CD1"))
            {$connect_x=substr($rec[$y],30,8);
             $connect_y=substr($rec[$y],38,8);
             $connect_z=substr($rec[$y],46,8);
             $vector3x=$connect_x-$x_r;
             $vector3y=$connect_y-$y_r;
             $vector3z=$connect_z-$z_r;
            }
         }
       $vector4x=$x_n[$k]-$x_r;
       $vector4y=$y_n[$k]-$y_r;
       $vector4z=$z_n[$k]-$z_r;
       $cos2=($vector3x*$vector4x+$vector3y*$vector4y+$vector3z*$vector4z)/(sqrt($vector3x**2+$vector3y**2+$vector3z**2)*sqrt($vector4x**2+$vector4y**2+$vector4z**2));
       $angle2=acos($cos2)*180/pi;
       $angle2=sprintf("%.4f",$angle2);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with $res :$dist_N A, with angle1 $angle1,angle2 $angle2\n";
      }
    elsif(($dist_N<3.4)&&((substr($rec[$z],17,3) eq "HOH")||(substr($rec[$z],17,3) eq "SOL")))
      {
       $vector2x=$x_r-$x_n[$k];
       $vector2y=$y_r-$y_n[$k];
       $vector2z=$z_r-$z_n[$k];
       $cos1=($vector1x*$vector2x+$vector1y*$vector2y+$vector1z*$vector2z)/(sqrt($vector1x**2+$vector1y**2+$vector1z**2)*sqrt($vector2x**2+$vector2y**2+$vector2z**2));
       $angle1=acos($cos1)*180/pi;
       $angle1=sprintf("%.4f",$angle1);
       printf "%9s:%2s:Nitrile %2d: with %10s:%9.4f A, with angle1 %9.4f, angle2 %9.4f\n",$ARGV[0],$chainID[$i],$k,$res,$dist_N,$angle1,$angle2;
       #print "$ARGV[0]:$chainID[$i]:Nitrile $k :with water $res :$dist_N A, with angle1 $angle1\n";
      }
   }
 }
}

