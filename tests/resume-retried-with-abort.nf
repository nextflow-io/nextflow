#!/usr/bin/env nextflow
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

process LONG_SLEEP {
	tag "long-sleep"

	input:
	val x

	output:
	stdout

	script:
	'''
	echo "LONG_SLEEP start"
	sleep 30
	echo "LONG_SLEEP done"
	'''
}

process SMALL_SLEEP_RETRY {
	tag "small-sleep-retry"

	errorStrategy 'retry'
	maxRetries 1

	input:
	val x

	output:
	stdout

	script:
	"""
	echo "SMALL_SLEEP_RETRY attempt: ${task.attempt}"
	sleep 5

	if [[ ${task.attempt} -eq 1 ]]; then
	  echo "Failing first attempt on purpose"
	  exit 1
	fi

	echo "Second attempt succeeded"
	"""
}

workflow {
	ch = Channel.of(1)
    LONG_SLEEP(ch)
	SMALL_SLEEP_RETRY(ch)
}
