#!/usr/bin/perl
use warnings;
use strict;

# Verification tool for checking vecgen-avx register assignments.
# Pass it the --dump-asm output from a single test.

# FIXME: this missed the case of a register whose reset at the end of
# a loop was missed.  Augment test.

my $logname = shift || "<stdin>";

my %inmem; # vr's in memory
my %physr; # vr to pr
my %virtr; # pr to vr

my %scratch; # vr to RBP offset

my $exitcode = 0;

my $line=0;
while(<>) {
    $line++;
    next if /^ /;
    chomp;
    if(/^FILL (.*) -> (.*)/) {
	my ($reg, $ymm) = ($1, $2);
	#print "FILL $ymm from $reg\n"; #DEBUG

	if(!$inmem{$reg} and $reg !~ /\/[A-Z]/) {
	    err("ERROR: filling unspilled register $reg to $ymm");
	}
	$physr{$reg} = $ymm;
	$virtr{$ymm} = $reg;

	#checkscratch("FILL", $reg);

    } elsif(/^SPILL (.*) <- (.*)/) {
	my ($reg, $ymm) = ($1, $2);
	#print "SPILL $ymm to $reg\n"; #DEBUG

	err("ERROR: spill of unfilled register $reg from $ymm")
	    if !defined $physr{$reg} or $physr{$reg} ne $ymm or $virtr{$ymm} ne $reg;
	$inmem{$reg}++;

	#checkscratch("SPILL", $reg);

    } elsif(/(.*) = ([A-Z]+)(.*)/) {
	my ($dst, $op, $args) = ($1, $2, $3);
	#print "$op dst: $dst\n"; #DEBUG

	for(my $argn=0; $args =~ /(R[0-9]+(?:\/[A-Z0-9]+)?)@(ymm[0-9]+)/g; $argn++) {
	    my ($reg, $ymm) = ($1, $2);
	    #print "  $reg @ $ymm (n $argn)\n"; #DEBUG

	    if((($op eq "LOAD" or $op eq "STORE")    and $argn != 0) or
	       (($op eq "CULL" or $op eq "CULL_FLD") and $argn != 1 )) {
		if(!$inmem{$reg}) {
		    if($reg !~ /\//) {
			# Only scratch (no slash) registers must be spilled
			err("ERROR: LOAD/STORE/CULL arg $reg not spilled");
		    }
		}
	    } else {
		# Everything else must be in registers
		if(!defined $virtr{$ymm} or $virtr{$ymm} ne $reg or
		   !defined $physr{$reg} or $physr{$reg} ne $ymm) {
		    err("ERROR: $op arg $reg not filled into $ymm");
		}
	    }
	}

	if($dst !~ /(R[0-9]+(?:\/[A-Z0-9]+)?)(?:@(ymm[0-9]+))?/) {
	    err("ERROR: unparsable dst: $dst");
	} else {
	    my ($reg, $ymm) = ($1, $2);
	    if($op eq "LOAD") {
		# LOADS deposit their results in memory
		delete $physr{$reg};
		$inmem{$dst}++;
	    } else {
		# All other destinations are registers, the in-memory
		# version is now stale
		delete $inmem{$dst};
		$virtr{$ymm} = $reg;
		$physr{$reg} = $ymm;
	    }
	}

    } elsif(/^LOOP/) {
	#print "LOOP\n"; #DEBUG

	# Mark register state for checking at POOL

    } elsif(/^POOLTEST (R[0-9]+(?:\/[A-Z0-9]+)?)@(ymm[0-9]+)/) {
	my ($reg, $ymm) = ($1, $2);
	#print "POOLTEST $reg in $ymm"; #DEBUG

	if(!defined $virtr{$ymm} or $virtr{$ymm} ne $reg or
	   !defined $physr{$reg} or $physr{$reg} ne $ymm) {
	    err("ERROR: POOLTEST arg $reg not filled into $ymm");
	}

    } elsif(/^POOL\(/) {
	#print "POOL\n"; #DEBUG

	# Check register state vs. the matching LOOP.
    }
}

# Note that this assumes the mapping of scratch location to vr is 1:1,
# that may change with future optimization.
#
# FIXME: disabled for now, not LOAD-cognizant (LOAD acts like a spill)
sub checkscratch {
    my ($mode, $r) = @_;
    return if $r =~ /\//; # not scratch
    my $s = $1 if <> =~ /(-0x[0-9a-f]+)\(%rbp\)/;
    $line++;
    if(!defined $s) {
	err("cannot find scratch offset for $r");
	return;
    }
    if(!defined $scratch{$r} && $mode eq "SPILL") {
	$scratch{$r} = $s; # first spill, not notable
    } elsif(!defined $scratch{$s}) {
	err("Attempt to $mode $r before first spill");
    } elsif($s != $scratch{$s}) {
	err("Multiple scratch offsets for $r: $s vs. $scratch{$s}");
    }
}

exit $exitcode;

sub err { print "$logname:$line: @_\n"; $exitcode = 1; }
