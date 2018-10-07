#!/usr/bin/perl -w
use strict;
use Math::CDF;
use Math::Round;

#
# The default input parameters 
# 
my $debug = 0;
my $noplot = 0;
my $gnuplot = "";
my $input_file = "protein.txt";
my $window_size = 7;
my $threshold = 0.225;
my $graph_file = "";
my $result_file = "results.$$";
my $data_file = "sliding.$$";
my $_prompt = "";

#
# Find out gnuplot 
#

if( exists $ENV{GNUPLOT_BIN}) { 
        $gnuplot = $ENV{GNUPLOT_BIN};
}
if( $gnuplot eq "" ) { 
        $gnuplot = `which gnuplot`;
        chomp($gnuplot);
}
if( $gnuplot eq "" ) { 
        print "AMPA requires gnuplot to be installed on your system.";
        exit 2;
}

#
# Parse the command 
# 

my $cl=join( " ", @ARGV);

if (($cl=~/-h/) ||($cl=~/-H/) ) {
        # print usage 
        print "Usage: AMPA.pl -in=<input fasta file> [other options]\n";
        print "\n";
        print "Available options:\n";
        print "-in      Fasta sequence input file\n";
        print "-t       Threshold value (default: 7)\n";
        print "-w       Window size value (default: 0.225)\n";
        print "-rf      File where store the program result\n";
        print "-df      File where store the produced plot data\n";
        print "-gf      File where store the generated plot\n";
        print "-noplot  Skip plot creation step\n";
        exit 1;
}

#
# If not command line is provided switch to interactive mode 
#

if ($#ARGV==-1 ) { 
        # the input file name
        print "Enter the name of the sequence file and then press Enter [$input_file]: ";
        $_prompt = <STDIN>;
        chomp($_prompt);
        if( $_prompt ne "" ) { $input_file=$_prompt };

        # the window size
        print "Enter the desired window size as an odd number and then press Enter [$window_size]: ";
        $_prompt = <STDIN>;
        chomp($_prompt);
        if( $_prompt ne "" ) { $window_size = $_prompt };
        
        # the threshold
        print "Enter the desired threshold and then press Enter [$threshold]: ";
        $_prompt = <STDIN>;
        chomp($_prompt);
        if( $_prompt ne "" ) { $threshold = $_prompt };
}


# option '-in': input file name
if ( ($cl=~/-in=\s*(\S+)/)) { 
        $input_file = $1;
}

# option '-w': window size
if ( ($cl =~ /-w=\s*(\S+)/)) { 
        $window_size = $1;
}

# option '-t': threshold value
if ( ($cl =~ /-t=\s*(\S+)/)) { 
        $threshold = $1;
}

# option '-gf': output graph file
if ( ($cl =~ /-gf=\s*(\S+)/)) { 
        $graph_file = $1;
}

# option '-rf': output result file
if ( ($cl =~ /-rf=\s*(\S+)/)) { 
        $result_file = $1;
}

# option '-df': output sliding file
if ( ($cl =~ /-df=\s*(\S+)/)) { 
        $data_file = $1;
}

# option '-debug': enable debug mode
if ( ($cl =~ /-debug/)) { 
        $debug = 1;
}

# option '-noplot': threshold value
if ( ($cl =~ /-noplot/)) { 
        $noplot = 1;
}

if( $debug != 0 ) { 
        print "Window size      : $window_size\n";
        print "Threshold value  : $threshold\n";
        print "Input file name  : $input_file\n";
        print "Graph file name  : $graph_file\n";
        print "Result file name : $result_file\n";
        print "Data file name   : $data_file\n";
        print "Gnuplot bin      : $gnuplot\n";
}

#
#  Main script start here 
# 

my(@names,@lengths,@sequences,$sequence_number,$protein);

open(PROTEIN, "$input_file") or die "Cannot open '$input_file': $!\n";  

#The program should read the user input sequence

$/ = ">";

for(my $c=1; $c<(($window_size + 1)/2); $c++) {
    push (@sequences, 'X');
    }

while (my $sequence_record = <PROTEIN>) {
    if ($sequence_record eq ">") {
        next;
    }

    my $sequence_name = "";
    if ($sequence_record =~ m/([^\n]+)/) {
        $sequence_name = $1;
    } else {
        $sequence_name = "No title was found!";
    }

    $sequence_record =~ s/[^\n]+//;
    push (@names, $sequence_name);
    $sequence_record =~ s/[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]//g;    
    push (@sequences, $sequence_record);
    push (@lengths, length($sequence_record));
    my $sequence_length = length($sequence_record);
}

for(my $c=1; $c<(($window_size + 1)/2); $c++) {
    push (@sequences, 'X');
}

my $sequences = join("",@sequences);

close(PROTEIN) or die "Cannot close the file: $!\n";

#The program calls the subroutines to analyze the sequence

&sliding($sequences);
if(!$noplot) {
	&gnuplot_sld();
} 

### Subroutines

# Defines the value at the window center (window_size+1)/2

sub sliding {
    open(KD,">$data_file") or die("Could not open '$data_file': $!\n");
    open(RS,">$result_file") or die("Could not open '$result_file': $!\n");
    my($center,$length);
    my $half = (($window_size + 1)/2);

# These are the antimicrobial propensity values for each amino acid

    my %sliding = (
                          'A',    0.307,  
                          'R',    0.106,  
                          'N',    0.240,  
                          'D',    0.479,  
                          'C',    0.165,  
                          'Q',    0.248,  
                          'E',    0.449,  
                          'G',    0.265,  
                          'H',    0.202,  
                          'I',    0.198,  
                          'L',    0.246,  
                          'K',    0.111,  
                          'M',    0.265,  
                          'F',    0.246,  
                          'P',    0.327,  
                          'S',    0.281,  
                          'T',    0.242,  
                          'W',    0.172,  
                          'Y',    0.185,  
                          'V',    0.200,
              'X',      $threshold  
                          );

    $protein = shift;

    print KD "# The window size is currently set at $window_size.\n";
    print KD "# Here are the AMSI values for this protein:\n";
    print KD "# $names[0]\n";
    print KD "# Pos\ APV\n";
    print KD "# ---\t-------------------\n";
        
        my $acc1=0;
        my $acc2=0;
        my $bacindex=0;
        my $amvalue=0;
    my $bacvalue=0;
        my $prob=0;

    for(my $i=0;$i<=(length($protein)-$window_size);$i++) {
        my $window = substr($protein, $i, $window_size);
        my $sum=0;
        for(my $j=0; $j<$window_size; $j++) {
            my $PV = 0;
            my $residue = substr($window, $j, 1);
            $PV = $sliding{$residue};
            $sum+=$PV;
        }
    
        $center = $i + $half;
    my $position = $center - 3;
        my $APV=$sum/$window_size;
        $APV = sprintf "%.3f", $APV;
        print KD "$position\t$APV\n";
        $amvalue=$amvalue+$APV;

                                if($acc1==3) {
                                                if($acc2>=10){
                                                        my $init=$position-$acc2-2;
                                                        my $last=$position-1;
                            my $bacvalue=(($bacvalue)/($last-$init+1));
                            $prob = &Math::CDF::pnorm(($bacvalue-0.2584)/0.02479);
                            $prob*=100;
                            my $rbacvalue = nearest(.001, $bacvalue);                            
                            my $rprob = nearest(1, $prob);
                                                        print RS "Antimicrobial stretch found in $init to $last. Propensity value $rbacvalue ($rprob %) \n";
                            $bacindex++;                                                        
                            $acc2=0;
                                                        $acc1=0;
                            $bacvalue=0;
                            $prob=0;
                                                        }
                                                if($acc2<10){
                                                $acc1=0;
                                                $acc2=0;
                        $bacvalue=0;
                                                }
                                             }
                                if($acc1!=3){
                        if($acc2==0) {
                            $acc1=0;
                            $bacvalue=0;
                        }                                                
                        if($APV<=$threshold) {
                                                $acc2++;
                        $bacvalue+=$APV
                                                }
                                                if($APV>$threshold){
                                                $acc1++;
                        $bacvalue+=$APV
                                                }
                                            }
                if($position==(length($protein)-(($window_size)-1)) && $acc2>=10) {
                        my $init=$position-$acc2-2;
                                                my $last=$position-1;
                        my $bacvalue=(($bacvalue)/($last-$init+1));
                        $bacindex++;
                        $prob = &Math::CDF::pnorm(($bacvalue-0.2584)/0.02479);
                        $prob*=100;
                        my $rbacvalue = nearest(.001, $bacvalue);                            
                        my $rprob = nearest(1, $prob);                                                
                        print RS "Antimicrobial stretch found in $init to $last. Propensity value $rbacvalue ($rprob %) \n";    
                        }


}


print RS "# This protein has $bacindex bactericidal stretches \n";
my $manvalue = nearest(.001, $amvalue/(length($protein)-6));
print RS "# This protein has a mean antimicrobial value of $manvalue \n";       
   
close(KD) or die("Could not close file: $!\n");
}


sub gnuplot_sld {
    my $cmdline = "|$gnuplot -persist";
   
    open (GP, $cmdline) or die "no gnuplot";
    # force buffer to flush after each write
    use FileHandle;
    GP->autoflush(1);

    print GP "set title \"Antimicrobial Profile\"\n";
    print GP "set xlabel \"Position\"\n";
    print GP "set ylabel \"Score\"\n";
    print GP "set grid\n";
    print GP "set data style lines\n";

        # if an outpout has been specified output to it 
    if( $graph_file ne "" ) { 
        print GP "set terminal png large enhanced  size 800 600\n";
        print GP "set output '$graph_file' \n";
    }
    
    print GP "plot \"$data_file\" u 1:2 t \"Antimicrobial profile \" w lines\n";
    close GP;
}

