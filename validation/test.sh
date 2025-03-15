#!/bin/bash 
set -e 

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_IGNORE_WARN_DSL2=true
export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}
# disable ansi log to make log more readable
export NXF_ANSI_LOG=false
export NXF_DISABLE_CHECK_LATEST=true

test_integration() {
    (
      cd "$1"
      sudo bash cleanup.sh
      cd checks
      bash run.sh
    )
}

test_e2e() {
    if [[ $TEST_SMOKE == true ]]; then
      echo Skipping tests since TEST_SMOKE flag is true
      exit 0
    fi

    #
    # Hello
    #
    git clone https://github.com/nextflow-io/hello
    (
      cd hello
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
}

#
# Integration tests
#
if [[ $TEST_MODE == 'test_integration' ]]; then
  test_integration ../tests/
  test_integration ../tests-v1/
  test_e2e
fi

#
# Integration tests (strict syntax)
#
if [[ $TEST_MODE == 'test_parser_v2' ]]; then
  export NXF_SYNTAX_PARSER=v2
  test_integration ../tests/
  test_e2e
fi

#
# Documentation tests
#
if [[ $TEST_MODE == 'test_docs' ]]; then

    (
      echo "Documentation tests"
      cd ../docs/snippets/
      bash test.sh
    )

fi

#
# AWS Batch tests
#
if [[ $TEST_MODE == 'test_aws' ]]; then
    if [ "$AWS_ACCESS_KEY_ID" ]; then
      echo "AWS batch tests"
      bash awsbatch.sh
    else
      echo "::warning file=$0,line=$LINENO::Missing AWS_ACCESS_KEY_ID variable -- Skipping AWS Batch tests"
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
      echo "::warning file=$0,line=$LINENO::Missing AZURE_BATCH_ACCOUNT_KEY variable -- Skipping Azure Batch tests"
    fi
fi

#
# Google Batch
#
if [[ $TEST_MODE == 'test_google' ]]; then
    if [ "$GOOGLE_SECRET" ]; then
      echo "Google Batch tests"
      bash google.sh
    else
      echo "::warning file=$0,line=$LINENO::Missing GOOGLE_SECRET variable -- Skipping Google Batch tests"
    fi
fi

#
# Wave
#
if [[ $TEST_MODE == 'test_wave' ]]; then
    if [ "$TOWER_ACCESS_TOKEN" ]; then
      echo "Wave tests"
      bash wave.sh
    else
      echo "::warning file=$0,line=$LINENO::Missing TOWER_ACCESS_TOKEN variable -- Skipping Wave tests"
    fi
fi
