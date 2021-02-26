#!/usr/bin/perl
open(LIST,"list");
open(NEW,"list_identify_CN");
@list=<LIST>;
@new=<NEW>;
foreach $line(@list)
  {$i=0;
   $pdbid=substr($line,0,4);
   foreach $line2(@new)
     {
      $newid=substr($line2,0,4);
      if ($pdbid eq $newid)
        {
         $i++;
        }
      }
    if($i==0)
      {
       print "$line";
      }
  }
