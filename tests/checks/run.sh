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
REPORT=$PWD/.report
WITH_DOCKER=${WITH_DOCKER:=''}
TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST:=false}

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
    if [[ `grep -c "$1" .IGNORE` != 0 ]]; then
        echo 'no'
    elif [[ ! $WITH_DOCKER && `grep -c "$1" .IGNORE-DOCKER` != 0 ]]; then
        echo 'no'
    elif [[ $TRAVIS_PULL_REQUEST != false && `grep -c "$1" .IGNORE-TRAVIS-PR` != 0 ]]; then 
        # https://docs.travis-ci.com/user/pull-requests/#Pull-Requests-and-Security-Restrictions
        echo 'no'    
    elif [[ -f .IGNORE-JAVA-$TEST_JDK && `grep -c "$1" .IGNORE-JAVA-$TEST_JDK` != 0 ]]; then
        echo 'no'
    else
        echo 'yes'
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
