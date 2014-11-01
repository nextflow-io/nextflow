#!/bin/bash
#
#  Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
#
#  This file is part of Nextflow.
#
#  Nextflow is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Nextflow is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.

#
# test single class with Gradle 
#
# Read more 
# http://www.gradle.org/docs/1.10/release-notes#executing-specific-tests-from-the-command-line
# http://stackoverflow.com/questions/18061774/run-single-integration-test-with-gradle
#

MODULE=${2:-''}

if [[ $# -eq 0 ]]; then
  echo 'usage: g-test <ClassToTest> [sub-project]'
  echo ''
  echo 'examples:'
  echo '  //select specific test method'
  echo '  gradle test --tests org.gradle.SomeTest.someFeature'
  echo ''
  echo '  //select specific test class'
  echo '  gradle test --tests org.gradle.SomeTest'
  echo ''
  echo '  //select all tests from package'
  echo '  gradle test --tests org.gradle.internal*'
  echo ''
  echo '  //select all ui test methods from integration tests by naming convention'
  echo '  gradle test --tests *IntegTest*ui*'
  echo ''
  echo '  //selecting tests from different test tasks'
  echo '  gradle test --tests *UiTest integTest --tests *WebTest*ui'
  echo ''
  exit 1
fi 

set -x
./gradlew -q $MODULE:test --tests $1
