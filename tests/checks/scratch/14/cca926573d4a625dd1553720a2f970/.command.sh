#!/usr/bin/env perl
    use strict;
    use warnings;

    my $count;
    my $range = 100;
    for ($count = 0; $count < 10; $count++) {
     	print rand($range) . ', ' . rand($range) . "
";
    }
