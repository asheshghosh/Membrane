#!/usr/bin/perl -w

use Math::Trig;
use Math::Round;
#use strict;
#require 'calcVecs.pl';

# Generate membrane surface image with protein chains
# Written by Ashesh Ghosh 07-07-2023

# Initialize the membrane parameters

my $kap=0.008*0.8;
my $sig=0.01*0.8;
$kap=0.005*0.8/1.7;
$sig=0.02*0.8/1.7;

my $l=250;
my $qmax=25;
my $pi=3.14159;
my $ngrid=200;
my $zcdel=2.5;
my $randc=1.2;
my $nphi=40;
my $nth=30;

# Define parameters for brush component 1

my $eps=800; # Rigidity of the brush polymers
$eps=500/10; # Rigidity of the brush polymers to be rendered
my $np=40; # Number of segments per polymer
my $l0p=2.5; # length of the segments
$l0p=1.2; # length of the segments
my $rho=0.025; # Density of polymers on membrane surface
my $frigid=0.80; #Density of rigid polymers
$rhor = $frigid*$rho;
$rhof = (1.-$frigid)*$rho;

# Variables in the calculation

my $ii;
my $jj;
my $kk;
my $qx;
my $qy;

# Calculate the surface

my $xs;
my $ys;
my $zs;

my $count=1;
$ii=1;
while ($ii <= $ngrid) {
  $jj=1;
  while ($jj <= $ngrid) {
    $xs[$count]=$l*($ii-1)/($ngrid-1)-$l/2;
    $ys[$count]=$l*($jj-1)/($ngrid-1)-$l/2;
    $zs[$count]=0;
    
    $count++;
    $jj++; 
  }
  $ii++;
}

my $ns=$count-1;

my $hq;
my $have;
my $qmag;

$qx=1;
while ($qx <= $qmax) {
  $qy=1;
  while ($qy <= $qmax) {
    $qmag=sqrt($qx**2+$qy**2);
    $have=1/($kap*$qmag**4+$sig*$qmag**2);
    $hq=$have*(rand()-0.5);
    $ii=1;
    $count=1;
    while ($ii <= $ngrid) {
      $jj=1;
      while ($jj <= $ngrid) {
	$zs[$count]=$zs[$count]+$hq*cos(2*$pi*($qx*$xs[$count]+$qy*$ys[$count])/$l);
	$count++;
	$jj++; 
      }
      $ii++;
    }
    $qy++;
  }
  $qx++;
}


# Calculate the brush layer

my $xb;
my $yb;
my $zb;
my $xb2;
my $yb2;
my $zb2;
my $nb=0;
my $nbtot1=0;
my $nbtot2=0;
my $theta;
my $phi;

my $t1x;
my $t1y;
my $t1z;
my $t2x;
my $t2y;
my $t2z;
my $t3x;
my $t3y;
my $t3z;

my $t1xp;
my $t1yp;
my $t1zp;
my $t2xp;
my $t2yp;
my $t2zp;
my $t3xp;
my $t3yp;
my $t3zp;


#first type of polymer

$count=1;
$countp=1;
$ii=1;
while ($ii <= $ngrid) {
    $jj=1;
    while ($jj <= $ngrid) {

        if (rand() <= $rhor) {
            cross:
            
            $nb++;
            $xb[$countp]=$xs[$count];
            $yb[$countp]=$ys[$count];
            $zb[$countp]=$zs[$count]+$l0p;
            $countp++;
            
            $t1x=1;
            $t1y=0;
            $t1z=0;
            
            $t2x=0;
            $t2y=1;
            $t2z=0;
            
            $t3x=0;
            $t3y=0;
            $t3z=1;
            
            $eps = 800/10;

            $kk=1;
            
            while ($kk <= ($np-1)) {
                
                $theta=acos(1/$eps*log(rand()*(exp($eps)-exp(-$eps))+exp(-$eps)));
                $phi=2*$pi*rand();
                
                $xb[$countp]=$xb[$countp-1]+$l0p*$t3x;
                $yb[$countp]=$yb[$countp-1]+$l0p*$t3y;
                $zb[$countp]=$zb[$countp-1]+$l0p*$t3z;
                if ($zb[$countp] <= ($zs[$count]+$l0p)) {
                    print "cross!\n";
                    $countp=$countp-$kk;
                    goto cross;
                }
                $countp++;
                
                
                $t1xp=cos($theta)*cos($phi)*$t1x+cos($theta)*sin($phi)*$t2x-sin($theta)*$t3x;
                $t1yp=cos($theta)*cos($phi)*$t1y+cos($theta)*sin($phi)*$t2y-sin($theta)*$t3y;
                $t1zp=cos($theta)*cos($phi)*$t1z+cos($theta)*sin($phi)*$t2z-sin($theta)*$t3z;
                    
                $t2xp=-sin($phi)*$t1x+cos($phi)*$t2x;
                $t2yp=-sin($phi)*$t1y+cos($phi)*$t2y;
                $t2zp=-sin($phi)*$t1z+cos($phi)*$t2z;
                
                $t3xp=sin($theta)*cos($phi)*$t1x+sin($theta)*sin($phi)*$t2x+cos($theta)*$t3x;
                $t3yp=sin($theta)*cos($phi)*$t1y+sin($theta)*sin($phi)*$t2y+cos($theta)*$t3y;
                $t3zp=sin($theta)*cos($phi)*$t1z+sin($theta)*sin($phi)*$t2z+cos($theta)*$t3z;
                
                $t1x=$t1xp;
                $t1y=$t1yp;
                $t1z=$t1zp;
                
                $t2x=$t2xp;
                $t2y=$t2yp;
                $t2z=$t2zp;
                
                $t3x=$t3xp;
                $t3y=$t3yp;
                $t3z=$t3zp;
                
                
                $kk++;
            }
        }
        $jj++;
        $count++;
    }
    $ii++;
}
$nbtot1=$countp-1;


#second type of polymer
$count=1;
$countp=1;
$ii=1;
while ($ii <= $ngrid) {
    $jj=1;
    while ($jj <= $ngrid) {

        if (rand() <= $rhof) {
            cross:
            
            $nb++;
            $xb2[$countp]=$xs[$count];
            $yb2[$countp]=$ys[$count];
            $zb2[$countp]=$zs[$count]+$l0p;
            $countp++;
            
            $t1x=1;
            $t1y=0;
            $t1z=0;
            
            $t2x=0;
            $t2y=1;
            $t2z=0;
            
            $t3x=0;
            $t3y=0;
            $t3z=1;
            
            $eps = 100/20;

            $kk=1;
            
            while ($kk <= ($np-1)) {
                
                $theta=acos(1/$eps*log(rand()*(exp($eps)-exp(-$eps))+exp(-$eps)));
                $phi=2*$pi*rand();
                
                $xb2[$countp]=$xb2[$countp-1]+$l0p*$t3x;
                $yb2[$countp]=$yb2[$countp-1]+$l0p*$t3y;
                $zb2[$countp]=$zb2[$countp-1]+$l0p*$t3z;
                if ($zb2[$countp] <= ($zs[$count]+$l0p)) {
                    print "cross!\n";
                    $countp=$countp-$kk;
                    goto cross;
                }
                $countp++;
                
                
                $t1xp=cos($theta)*cos($phi)*$t1x+cos($theta)*sin($phi)*$t2x-sin($theta)*$t3x;
                $t1yp=cos($theta)*cos($phi)*$t1y+cos($theta)*sin($phi)*$t2y-sin($theta)*$t3y;
                $t1zp=cos($theta)*cos($phi)*$t1z+cos($theta)*sin($phi)*$t2z-sin($theta)*$t3z;
                    
                $t2xp=-sin($phi)*$t1x+cos($phi)*$t2x;
                $t2yp=-sin($phi)*$t1y+cos($phi)*$t2y;
                $t2zp=-sin($phi)*$t1z+cos($phi)*$t2z;
                
                $t3xp=sin($theta)*cos($phi)*$t1x+sin($theta)*sin($phi)*$t2x+cos($theta)*$t3x;
                $t3yp=sin($theta)*cos($phi)*$t1y+sin($theta)*sin($phi)*$t2y+cos($theta)*$t3y;
                $t3zp=sin($theta)*cos($phi)*$t1z+sin($theta)*sin($phi)*$t2z+cos($theta)*$t3z;
                
                $t1x=$t1xp;
                $t1y=$t1yp;
                $t1z=$t1zp;
                
                $t2x=$t2xp;
                $t2y=$t2yp;
                $t2z=$t2zp;
                
                $t3x=$t3xp;
                $t3y=$t3yp;
                $t3z=$t3zp;
                
                
                $kk++;
            }
        }
        $jj++;
        $count++;
    }
    $ii++;
}
$nbtot2=$countp-1;



# Assemble single PDB file

$fileout=">membrane.pdb";
open(PDB, $fileout) || die('cannot open file:'. $!);

my $atomname1 = "A1";           # Chain atom type
my $atomname2 = "A2";           # Ribbon atom type
my $atomname3 = "A3";           # Extra atom type
my $atomname4 = "A4";           # Extra atom type
my $atomname5 = "A5";           # Extra atom type
my $atomname6 = "A6";           # Extra atom type
my $atomname7 = "A7";           # Extra atom type
my $atomname8 = "A8";           # Extra atom type
my $resname = "SSN";           # Type of residue (UNKnown/Single Stranded Nucleotide)
my $chain = "A";               # Chain identifier
my $resnum = "1";
my $numresidues = 1;
my $descrip = "Pseudo atom representation of DNA";
my $chemicalname = "Body and ribbon spatial coordinates";

# Het Header info
printf PDB "HET    %3s  %1s%4d   %5d     %-38s\n",$resname,$chain,$resnum, $numresidues,$descrip ;
printf PDB "HETNAM     %3s %-50s\n",$resname, $chemicalname;
printf PDB "FORMUL  1   %3s    C20 N20 P21\n",$resname;

$count=1;
$ii=1;
while ($ii <= $ns) {
  printf PDB "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",$count,$atomname1,$resname,1,$xs[$ii],$ys[$ii],$zs[$ii],1,1;
  $count++;
  $ii++;
}

$ii=1;
while ($ii <= $nbtot1) {
    printf PDB "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",$count,$atomname2,$resname,1,$xb[$ii],$yb[$ii],$zb[$ii],1,1;
    $count++;
    $ii++;
}

$ii=1;
while ($ii <= $nbtot2) {
    printf PDB "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",$count,$atomname3,$resname,1,$xb2[$ii],$yb2[$ii],$zb2[$ii],1,1;
    $count++;
    $ii++;
}

# Connect the grid points together

$ii=1;
while ($ii <= $ngrid) {
  $jj=1;
  while ($jj <= $ngrid) {
    if ($ii==1) {
      if ($jj==1) {
	printf PDB "CONECT%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1)+$jj+1;
      } elsif ($jj==$ngrid) {
	printf PDB "CONECT%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1)+$jj-1;
      } else {
	printf PDB "CONECT%5d%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1)+$jj+1,$ngrid*($ii-1)+$jj-1;
      }
    } elsif ($ii==$ngrid) {
      if ($jj==1) {
	printf PDB "CONECT%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj+1;
      } elsif ($jj==$ngrid) {
	printf PDB "CONECT%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj-1;
      } else {
	printf PDB "CONECT%5d%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj+1,$ngrid*($ii-1)+$jj-1;
      }
    } else {
      if ($jj==1) {
	printf PDB "CONECT%5d%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj+1;
      } elsif ($jj==$ngrid) {
	printf PDB "CONECT%5d%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj-1;
      } else {
	printf PDB "CONECT%5d%5d%5d%5d%5d\n",$ngrid*($ii-1)+$jj,$ngrid*($ii-1+1)+$jj,$ngrid*($ii-1-1)+$jj,$ngrid*($ii-1)+$jj+1,$ngrid*($ii-1)+$jj-1;
      }
    }
    $jj++;
  }
  $ii++;
}

$ii=1;
$ind0=$ns;
$count=1;
while ($ii <= $nbtot1) {
    if ($count==1) {
        printf PDB "CONECT%5d%5d\n",$ind0+$ii,$ind0+$ii+1;
        $count++;
    } elsif ($count==$np) {
        printf PDB "CONECT%5d%5d\n",$ind0+$ii,$ind0+$ii-1;
        $count=1;
    } else {
        printf PDB "CONECT%5d%5d%5d\n",$ind0+$ii,$ind0+$ii-1,$ind0+$ii+1;
        $count++;
    }
    
    $ii++;
}
while ($ii <= $nbtot2) {
    if ($count==1) {
        printf PDB "CONECT%5d%5d\n",$ind0+$ii,$ind0+$ii+1;
        $count++;
    } elsif ($count==$np) {
        printf PDB "CONECT%5d%5d\n",$ind0+$ii,$ind0+$ii-1;
        $count=1;
    } else {
        printf PDB "CONECT%5d%5d%5d\n",$ind0+$ii,$ind0+$ii-1,$ind0+$ii+1;
        $count++;
    }
    
    $ii++;
}

printf PDB "END";

close(PDB);
