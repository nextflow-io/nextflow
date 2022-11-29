#!/bin/bash 
set -e 

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_IGNORE_WARN_DSL2=true
export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}
# disable ansi log to make log more readable
export NXF_ANSI_LOG=false

#
# Integration test
#
if [[ $TEST_MODE == 'test_integration' ]]; then

    (
      cd ../tests/
      sudo bash cleanup.sh
      cd checks
      bash run.sh
    )

    if [[ $TEST_SMOKE == true ]]; then
      echo Skipping tests since TEST_SMOKE flag is true
      exit 0
    fi

    #
    # Hello
    #
    git clone https://github.com/nextflow-io/hello
    (
      cd hello;
      $NXF_CMD run .
      $NXF_CMD run . -resume
    )

    #
    # RNASEQ-NF
    #
    echo nextflow-io/rnaseq-nf
    [[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
    $NXF_CMD run nextflow-io/rnaseq-nf -with-docker $OPTS
    $NXF_CMD run nextflow-io/rnaseq-nf -with-docker $OPTS -resume

    exit 0
fi


if [[ $GITHUB_EVENT_NAME == pull_request ]]; then
  echo "Skipping cloud integration tests on PR event"
  exit 0
fi


#
# AWS Batch tests
#
if [[ $TEST_MODE == 'test_aws' ]]; then
    if [ "$AWS_ACCESS_KEY_ID" ]; then
      echo "AWS batch tests"
      bash awsbatch.sh
    else
      echo "Missing AWS_ACCESS_KEY_ID variable -- Skipping AWS Batch tests"
    fi
fi

#
# Azure Batch tests
#
if [[ $TEST_MODE == 'test_azure' ]]; then
    if [ "$AZURE_BATCH_ACCOUNT_KEY" ]; then
      echo "Azure batch tests"
      bash azure.sh
    else
      echo "Missing AZURE_BATCH_ACCOUNT_KEY variable -- Skipping Azure Batch tests"
    fi
fi

#
# Google Life Sciences
#
if [[ $TEST_MODE == 'test_google' ]]; then
    if [ "$GOOGLE_SECRET" ]; then
      echo "Google LS tests"
      bash google.sh
    else
      echo "Missing GOOGLE_SECRET variable -- Skipping Google LS tests"
    fi
fi
