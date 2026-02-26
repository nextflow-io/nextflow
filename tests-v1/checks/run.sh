#!/bin/bash 
set -u 
trap "exit" INT

NORMAL=$(tput sgr0)
GREEN=$(tput setaf 2; tput bold)
YELLOW=$(tput setaf 3)
RED=$(tput setaf 1)

function echo_red() {
    echo -e "$RED$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

function echo_yellow() {
    echo -e "$YELLOW$*$NORMAL"
}


#
# Some vars 
#
NXF_CMD=${NXF_CMD:-nextflow}
NXF_SYNTAX_PARSER=v1
REPORT=$PWD/.report
WITH_DOCKER=${WITH_DOCKER:=''}

#
# Clean scratch dir 
#
export NXF_WORK=$PWD/scratch
rm -rf $NXF_WORK


function run_checks() {
  NXF_SCRIPT="../../$1"
  NXF_RUN="$NXF_CMD -q run $NXF_SCRIPT"
  export NXF_SCRIPT
  export NXF_CMD
  export NXF_RUN
  set +e   
            
  if [ -f $1/.checks ]; then 
     cd $basename; 
     rm -rf *
     bash -ex .checks &> checks.out
  else
     mkdir -p $1
     cd $1
     $NXF_RUN > checks.out
  fi
  
  ret=$?
  set -e
  if [[ $ret != 0 ]]; then 
    echo "~ Test '$1' run failed" >> $REPORT
    # dump error output
    [[ -s checks.out ]] && cat checks.out | sed 's/^/   /'>> $REPORT  
    echo '' >> $REPORT 
    # dump nextflow log file  
    [[ -f .nextflow.log ]] && cat .nextflow.log >> $REPORT
    echo '' >> $REPORT   
    exit 1
  fi  
}


rm -rf $REPORT

list=${1:-'../*.nf'}

function can_run() {
    if [[ $(grep -c "$1" .PARSER-V1) != 0 ]]; then
        echo 'yes'
    else
        echo 'no'
    fi
}

for x in $list; do
  basename=$(basename $x)
  if [[ $(can_run $basename) == 'yes' ]]; then
    echo "> Running test: $basename"
    ( set -e; 
      run_checks $basename 
    )
  else
    echo "- Ignoring test: $basename"
  fi  
done

if [[ -s $REPORT ]]; then
  echo -e "$RED"
  cat $REPORT
  echo -e "$NORMAL"
  exit 1
fi 
