#!/usr/bin/perl

#
# (c) 2009-2010 andre severo pereira gomes (andre.gomes@univ-lille1.fr)_
#
# summary_fscc.pl will use the thresholds used to define the Pm/Pi partition in order to
# verify the number of determinants (and their weights on the electronic state)
#   a) "in Pm" (=those with orbital energies within  [EHMAX,EPMIN])
#   b) "in Pi" (=those with orbital energies outside [EHMAX,EPMIN]) 
#   c) "mixed" (=those where one of the orbitals is in [EHMAX,EPMIN] while the other is not)
#
# the output can be redirected to a text file and then imported in a spreadsheet program and further processed, 
# e.g. have the energies sorted (now they are printed in the order they appear in the output) 
# and/or converted to cm-1 or eV (they are in atomic units now, and are relative to the reference state
# for a given sector - now only sector 11 has been used, but extending it is simple) 
#

$pm_lower      = -0.80;
$pm_higher     =  0.19;
$read_sector   =  0;
$totalstates   =  0;
$relevance_thr =  0.05;

while (<>) 
{
   if ( /CCIH EHMIN=([-.0-9e+]+), EHMAX=([-.0-9e+]+), EPMIN=([-.0-9e+]+), EPMAX=([-.0-9e+ ]+)/ )
   {
      $pm_lower = $2;
      $pm_higher = $3;
      $p_lower  = $1;
      $p_higher = $4;
      print "# Pm space defined via orbital energies between ",$pm_lower," and ",$pm_higher,"\n";
      print "# P  space defined via orbital energies between ",$p_lower," and ",$p_higher,"\n";
   }

   if ( /sector (\d+)/ )
   {
      if ( $1 == "11" )
      {
         $read_sector = 1;
      }
      else
      {
         $read_sector = 0;
      }
   }
 
   if ( $read_sector )
   {
      if ( /Irrep\s+([-+ a-zA-Z0-9]+)\s+State\s+(\d+)\s+([-.0-9]+)\s+([-.0-9]+)/ )
      {
#        print "BLA!\n";

         $totalstates++;

         $irrep[$totalstates]  = $1;
         $abs_en[$totalstates] = $3;
         $rel_en[$totalstates] = $3;

         $weightinpm[$totalstates] = 0; 
         $dets_in_pm[$totalstates] = 0;
         $dets_signi[$totalstates] = 0;

         $weightinmixedpipm[$totalstates] = 0;
         $dets_in_mixedpipm[$totalstates] = 0;
         $dets_signimixedpi[$totalstates] = 0; 
         
      }

#     print "read line : ",$_;
#     if ( /^\s+([-.0-9]+)\s+([ +-a-zA-Z0-9]+)\s*\#\s+\d+\s*\(\s*([-.0-9]+)\)\s*\-\>\s*([ +-a-zA-Z0-9]+)\s*\#\s+\d+\s*\(\s*([-.0-9]+)\)/ )
      if ( /^\s+([-.0-9]+)\s+([ +-a-zA-Z0-9]+).+\d+\s*\(\s*([-.0-9]+)\).+\-\>\s+([ +-a-zA-Z0-9]+).+\d+\s*\(\s*([-.0-9]+)\)/ )
      {
#        print $_;

         $coef_det = $1;
         $from_orb = $3;
         $to_orb   = $5;

         if ( ($from_orb < $pm_higher) && 
              ($from_orb > $pm_lower)  &&
              ($to_orb   < $pm_higher) && 
              ($to_orb   > $pm_lower)     ) 
         {
#            print "\n",$_;
#            print "     => det certainly in pm! coef ",$coef_det,"\n\n";
             $wt = $coef_det * $coef_det;
             $weightinpm[$totalstates] += $wt;
             $dets_in_pm[$totalstates]++;

             $dets_signi[$totalstates]++ if ( $wt > $relevance_thr );
             
         }

         if (  ($from_orb < $pm_higher) &&
               ($to_orb   > $pm_lower)  &&
               ( 
                 (($from_orb < $pm_lower)&&($to_orb   < $pm_higher)) || 
                 (($from_orb > $pm_lower)&&($to_orb   > $pm_higher)) 
               ) 
            ) 
         {
#            print "\n",$_;
#            print "     => det partly in pm! coef ",$coef_det,"\n\n";
             $wt = $coef_det * $coef_det;
             $weightinmixedpipm[$totalstates] += $wt;
             $dets_in_mixedpipm[$totalstates]++;

             $dets_signimixedpi[$totalstates]++ if ( $wt > $relevance_thr );
             
         }
   
      }
   }
}

print "\n\n Summary \n\n";

for ($i = 1; $i <= $totalstates; $i++)
{
   $irrep[$i]  =~ s/ //g;
   printf("energy   %12.6f    label   %8s    percentage_pm  %12.6f    dets_pm   %4d    significant_dets   %4d    percentage_mixed_pi_pm  %12.6f    dets_mixed_pi_pm   %4d    significant_dets_mixed_pi_pm  %4d\n",
          $rel_en[$i],$irrep[$i],$weightinpm[$i],$dets_in_pm[$i],$dets_signi[$i],$weightinmixedpipm[$i],$dets_in_mixedpipm[$i],$dets_signimixedpi[$i]);
}
