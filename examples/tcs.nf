#!/usr/bin/env nextflow

// define the input parameters
NUM_SET=1           // 100
action='genAln'     // genAln, genFilter, genTree, genML
//method=
RUN_MODE = 'debug'

METHODS="GS GR WR TG TS"
LONG_METHODS="$METHODS OA"

DB_FOLDER='./data'
DATA_FOLDER='./results'
LOG_FOLDER='./logs'


if( RUN_MODE == "debug" ) {
  tips=["tips32"]           // tips number
  divers=["asymmetric_0.5"]
  lens=["0400"]             // aln len
  alns=["MA","CL","PC"]     // MAFFT, ClustalW2, ProbCons
}
else if( RUN_MODE == "complete"  ) {
  tips=["tips32","tips64"]
  divers=["asymmetric_0.5", "asymmetric_1.0","asymmetric_2.0","symmetric_0.5","symmetric_1.0","symmetric_2.0"]
  lens=["0400","0800","1200","1600","3200"]
  alns=["MA","CL","PC"] // MAFFT, ClustalW2, ProbCons

} else {
  tips=["tips32","tips64"]
  divers=["asymmetric_0.5", "asymmetric_1.0", "asymmetric_2.0"]
  lens=["0400", "0800", "1200"]
  alns=["MA"]
}


genAlnInputData = channel()


for( tip in tips )
  for( diver in divers )
    for( len in lens )
	    for( aln in alns ) {

	        switch( action ) {
	        case "genAln":
            	    println "qsub -N $action -t 1-$NUM_SET ./wrapper-4-genAln $tip $diver $len $aln"
            	    genAlnInputData << [ tip, diver, len, aln, 1 ]
            	    break

            case "genFilter":
            	    //LOG="${LOG_FOLDER}/${action}_method-${method}_aln-${aln}-\$JOB_ID.o"
            	    println '$qsub_cmd -N $action -o $LOG -t 1-$NUM_SET ./wrapper-4-genFilter $tip $diver $len $aln $method'
                    break

            case "genTree":
                    //LOG="${LOG_FOLDER}/${action}_method-${method}_aln-${aln}-\$JOB_ID.o"
            	    println '$qsub_cmd -N $action -o $LOG -t 1-$NUM_SET ./wrapper-4-genTree $tip $diver $len $aln $method'
            	    break

            case "genML":
                    //LOG="${LOG_FOLDER}/${action}_method-${method}_aln-${aln}-\$JOB_ID.o"
            	    println '$qsub_cmd -N $action -o $LOG -t 1-$NUM_SET ./wrapper-4-genML $tip $diver $len $aln $method'
            	    break

            default:
                println "Invalid action: '$action'"

	        }


	    }

/*
task('genAln') {

    input genAlnInputData

    def reformat_cmd="t_coffee -other_pg seq_reformat"
    def mafft_cmd="mafft --quiet --nj"
    def clustal_cmd="clustalw"
    def probcons_cmd="probcons"

    def (tip_p, var_p, len, aln, set_p) = genAlnInputData

    def input_f="${DB_FOLDER}/01.Data/$tip_p/$var_p/$set_p/seqs/$set_p.$len.fa"
    def out_p="${DATA_FOLDER}/$tip_p/$var_p/len$len/set$set_p/OA/$aln"
    def fa_f="MSA.fa"
    def phylip_f="MSA.phylip"

    """
    if [ -e $input_f ]
    then
      echo "process $input_f"

      # Skip if already exists
      if [ -e $out_p/$phylip_f ]; then echo "SKIP genAln: $out_p/$phylip_f exists"; exit 0; fi

       cp $input_f tmp.fasta
       case $aln in
       "MA")
         $mafft_cmd tmp.fasta > tmp.aln 2> /dev/null
         ;;
       "CL")
         $clustal_cmd -infile=tmp.fasta -outfile=tmp.aln
         ;;
       "PC")
         $probcons_cmd tmp.fasta > tmp.aln 2> /dev/null
         ;;
       esac
       $reformat_cmd -in tmp.aln -output phylip > $phylip_f
       $reformat_cmd -in tmp.aln -output fasta > $fa_f

       # move result to target folder
       [ ! -e $out_p ] && mkdir -p $out_p
       mv $fa_f $out_p
       mv $phylip_f $out_p

       rm tmp.*
    else
      echo "ERROR: $input_f not exists"
    fi
    """

}
   */