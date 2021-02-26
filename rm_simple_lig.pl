#!/usr/bin/perl
open(LIG,"list_lig");
open(SIM,"list_simple_lig");
open(NEW,">list_lig_new");
@lig=<LIG>;
@sim=<SIM>;
foreach $line(@lig)
  {$i=0;
   $id=substr($line,0,3);
   #print "$line";
   foreach $line2(@sim)
     {$id_2=substr($line2,0,3);
      if($id eq $id_2)
        {
        # print "It is $id_2\n";
         $i++;
        }
      }
   if($i==0)
     {
      print NEW ("$line");
     }
   }
