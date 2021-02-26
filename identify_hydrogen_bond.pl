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
 if(@{$chainID[$i]}/2>=2)
 {
 #print "The number of CN in $chainID[$i] is $cn_num\n";
 $x_n=substr(${$chainID[$i]}[0],30,8);
 $y_n=substr(${$chainID[$i]}[0],38,8);
 $z_n=substr(${$chainID[$i]}[0],46,8);
 $x_c=substr(${$chainID[$i]}[1],30,8);
 $y_c=substr(${$chainID[$i]}[1],38,8);
 $z_c=substr(${$chainID[$i]}[1],46,8);
 $x_n2=substr(${$chainID[$i]}[2],30,8);
 $y_n2=substr(${$chainID[$i]}[2],38,8);
 $z_n2=substr(${$chainID[$i]}[2],46,8);
 $x_c2=substr(${$chainID[$i]}[3],30,8);
 $y_c2=substr(${$chainID[$i]}[3],38,8);
 $z_c2=substr(${$chainID[$i]}[3],46,8);
 foreach $line(@rec)
   {$res=substr($line,17,9);
    $x_r=substr($line,30,8);
    $y_r=substr($line,38,8);
    $z_r=substr($line,46,8);
    $dist_N=sprintf("%.4f",sqrt(($x_r-$x_n)**2+($y_r-$y_n)**2+($z_r-$z_n)**2));
    $dist_C=sprintf("%.4f",sqrt(($x_r-$x_c)**2+($y_r-$y_c)**2+($z_r-$z_c)**2));
    $dist_N2=sprintf("%.4f",sqrt(($x_r-$x_n2)**2+($y_r-$y_n2)**2+($z_r-$z_n2)**2));
    $dist_C2=sprintf("%.4f",sqrt(($x_r-$x_c2)**2+($y_r-$y_c2)**2+($z_r-$z_c2)**2));
    @tmp=sort{$a<=>$b}$dist_N,$dist_N2;
    
    if(($tmp[0]<3.2)&&(substr($line,13,3) eq "OH "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&(substr($line,13,3) eq "OG "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&(substr($line,13,3) eq "OG1"))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&(substr($line,13,3) eq "N  "))
      {print "$ARGV[0]:$chainID[$i]:with $res backbone NH:$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&((substr($line,13,3) eq "NH1")||(substr($line,13,3) eq "NH2")||(substr($line,13,3) eq "NE ")))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&(substr($line,13,3) eq "NZ "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&((substr($line,13,3) eq "ND2")||(substr($line,13,3) eq "NE2")))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&((substr($line,13,3) eq "NE2")||(substr($line,13,3) eq "ND1"))&&((substr($line,17,3) eq "HIS")||(substr($line,17,3) eq "HIE")||(substr($line,17,3) eq "HID")))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&(substr($line,13,3) eq "NE1"))
      {print "$ARGV[0]:$chainID[$i]:with $res :$tmp[0] A\n";
      }
    elsif(($tmp[0]<3.2)&&((substr($line,17,3) eq "HOH")||(substr($line,17,3) eq "SOL")))
      {print "$ARGV[0]:$chainID[$i]:with water $res :$tmp[0] A\n";
      }
   }
 }
 if(@{$chainID[$i]}/2==1)
 {
 $x_n=substr(${$chainID[$i]}[0],30,8);
 $y_n=substr(${$chainID[$i]}[0],38,8);
 $z_n=substr(${$chainID[$i]}[0],46,8);
 $x_c=substr(${$chainID[$i]}[1],30,8);
 $y_c=substr(${$chainID[$i]}[1],38,8);
 $z_c=substr(${$chainID[$i]}[1],46,8);
 foreach $line(@rec)
   {$res=substr($line,17,9);
    $x_r=substr($line,30,8);
    $y_r=substr($line,38,8);
    $z_r=substr($line,46,8);
    $dist_N=sprintf("%.4f",sqrt(($x_r-$x_n)**2+($y_r-$y_n)**2+($z_r-$z_n)**2));
    $dist_C=sprintf("%.4f",sqrt(($x_r-$x_c)**2+($y_r-$y_c)**2+($z_r-$z_c)**2));
    if(($dist_N<3.2)&&(substr($line,13,3) eq "OH "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&(substr($line,13,3) eq "OG "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&(substr($line,13,3) eq "OG1"))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&(substr($line,13,3) eq "N  "))
      {print "$ARGV[0]:$chainID[$i]:with $res backbone NH:$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&((substr($line,13,3) eq "NH1")||(substr($line,13,3) eq "NH2")||(substr($line,13,3) eq "NE ")))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&(substr($line,13,3) eq "NZ "))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&((substr($line,13,3) eq "ND2")||(substr($line,13,3) eq "NE2")))
      {print "$ARGV[0]:$chainID[$i]:with $res : $dist_N A\n";
      }
    elsif(($dist_N<3.2)&&((substr($line,13,3) eq "NE2")||(substr($line,13,3) eq "ND1"))&&((substr($line,17,3) eq "HIS")||(substr($line,17,3) eq "HIE")||(substr($line,17,3) eq "HID")))
      {print "$ARGV[0]:$chainID[$i]:with $res : $dist_N A\n";
      }
    elsif(($dist_N<3.2)&&(substr($line,13,3) eq "NE1"))
      {print "$ARGV[0]:$chainID[$i]:with $res :$dist_N A\n";
      }
    elsif(($dist_N<3.2)&&((substr($line,17,3) eq "HOH")||(substr($line,17,3) eq "SOL")))
      {print "$ARGV[0]:$chainID[$i]:with water $res :$dist_N A\n";
      }
   }
 }
}

