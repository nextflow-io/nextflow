#!/usr/bin/env nextflow

params.range = 100

randNums = channel()

/*
 * A trivial Perl script producing a list of numbers pair
 */
task {
    output '-': randNums
    shell 'perl'

    """
    use strict;
    use warnings;

    my $count;
    my $range = ${params.range};
    for ($count = 0; $count < 10; $count++) {
     	print rand($range) . ', ' . rand($range) . "\\n";
    }

    """

}


/*
 * A Python script task which parses the output of the previous script
 */
task {
    input '-': randNums
    echo true
    shell 'python'

    """
    import sys

    x = 0
    y = 0
    lines = 0
    for line in sys.stdin:
        items = line.strip().split(",")
        x = x+ float(items[0])
        y = y+ float(items[1])
        lines = lines+1

    print "avg: %s - %s" % ( x/lines, y/lines )

    """


}

