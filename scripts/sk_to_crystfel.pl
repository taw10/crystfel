#!/usr/bin/perl -w

use strict;

open(ROT, "< rotations.lst");
open(COORD, "< coords.dat");
open(OUT, "> cspad-feb2011.geom");

my $quadrant = 0;
printf(OUT "n_panels = 64\n\n");
my $x = 0;
my $y = 0;
my $p = 0;

my $global_x = 0.0;
my $global_y = 0.0;

while ( my $coord = <COORD> ) {

	my $rot = <ROT>;
	my $cx = 0.0;
	my $cy = 0.0;
	my $minx = $x*194;
	my $miny = $y*185;
	my $sc1;
	my $sc2;

	chomp $coord;
	chomp $rot;

	if ( $coord =~ /\(([0-9\.\-]+),\s([0-9\.\-]+)\)/ ) {
		$sc1 = $1;
		$sc2 = $2;
		printf("%f %f\n", $sc1, $sc2);

	} else {
		printf("!!!\n");
	}

	my $a = 388.0; # Fast scan
	my $b = 185.0; # Slow scan

	# FIXME: Remove the $b for versions of Stephan's code after 20th Feb
	$sc1 -= (1308.696-$b);
	$sc2 -= (980.3862-$b);

	my $sx = -$sc2;
	my $sy = -$sc1;

	printf(OUT "; Quadrant %i, asic %i\n", $quadrant, ($x%2)+2*$y);
	printf(OUT "%i/min_fs = %i\n", $p, $minx);
	printf(OUT "%i/min_ss = %i\n", $p, $miny);
	printf(OUT "%i/max_fs = %i\n", $p, ($x+1)*194-1);
	printf(OUT "%i/max_ss = %i\n", $p, ($y+1)*185-1);
	printf(OUT "%i/badrow_direction = -\n", $p);
	printf(OUT "%i/res = 9090.91\n", $p);
	printf(OUT "%i/peak_sep = 6.0\n", $p);
	printf(OUT "%i/clen = 77.0e-3\n", $p);
	if ( $rot == "0" ) {
		printf(OUT "%i/fs = -x\n", $p);
		printf(OUT "%i/ss = -y\n", $p);
		$cx = $sx;
		$cy = $sy;
	} elsif ( $rot == "90" ) {
		printf(OUT "%i/fs = +y\n", $p);
		printf(OUT "%i/ss = -x\n", $p);
		$cx = $sx;
		$cy = $sy - $a - 5.0;
	} elsif ( $rot == "180" ) {
		printf(OUT "%i/fs = +x\n", $p);
		printf(OUT "%i/ss = +y\n", $p);
		$cx = $sx - $a - 5.0;
		$cy = $sy - $b;
	} elsif ( $rot == "270" ) {
		printf(OUT "%i/fs = -y\n", $p);
		printf(OUT "%i/ss = +x\n", $p);
		$cx = $sx - $b;
		$cy = $sy;
	}
	printf(OUT "%i/corner_x = %5.2f\n", $p, $cx+$global_x);
	printf(OUT "%i/corner_y = %5.2f\n", $p, $cy+$global_y);
	printf(OUT "%i/no_index = 0\n", $p);
	printf(OUT "\n");
	#if ( ($y < 20) &&  ($x < 1) ) {
		printf(STDERR "%f %f %f %f\n", $cx, $cy, $sc1, $sc2);
	#}
	$x++;
	$p++;

	printf(OUT "; Quadrant %i, asic %i\n", $quadrant, ($x%2)+2*$y);
	printf(OUT "%i/min_fs = %i\n", $p, $x*194);
	printf(OUT "%i/min_ss = %i\n", $p, $y*185);
	printf(OUT "%i/max_fs = %i\n", $p, ($x+1)*194-1);
	printf(OUT "%i/max_ss = %i\n", $p, ($y+1)*185-1);
	printf(OUT "%i/badrow_direction = -\n", $p);
	printf(OUT "%i/res = 9090.91\n", $p);
	printf(OUT "%i/peak_sep = 6.0\n", $p);
	printf(OUT "%i/clen = 77.0e-3\n", $p);
	if ( $rot == "0" ) {
		printf(OUT "%i/fs = -x\n", $p);
		printf(OUT "%i/ss = -y\n", $p);
		$cx = $sx - $a/2.0 - 5.0;
		$cy = $sy;
	} elsif ( $rot == "90" ) {
		printf(OUT "%i/fs = +y\n", $p);
		printf(OUT "%i/ss = -x\n", $p);
		$cx = $sx;
		$cy = $sy - $a/2.0;
	} elsif ( $rot == "180" ) {
		printf(OUT "%i/fs = +x\n", $p);
		printf(OUT "%i/ss = +y\n", $p);
		$cx = $sx - $a/2.0;
		$cy = $sy - $b;
	} elsif ( $rot == "270" ) {
		printf(OUT "%i/fs = -y\n", $p);
		printf(OUT "%i/ss = +x\n", $p);
		$cx = $sx - $b;
		$cy = $sy - $a/2.0 - 5.0;
	}
	printf(OUT "%i/corner_x = %5.2f\n", $p, $cx+$global_x);
	printf(OUT "%i/corner_y = %5.2f\n", $p, $cy+$global_y);
	printf(OUT "%i/no_index = 0\n", $p);
	printf(OUT "\n");
	#if ( ($y < 20) &&  ($x < 1) ) {
		printf(STDERR "%f %f %f %f\n", $cx, $cy, $sc1, $sc2);
	#}
	$x++;
	$p++;

	if ( $x == 8 ) {
		$x = 0;
		$y++;
	}

	$quadrant++;
	if ( $quadrant == 4 ) {
		$quadrant = 0;
	}
}
