#!/bin/bash
#
#  Copyright 2013-2019, Centre for Genomic Regulation (CRG)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

#
# This script acts as a pass-through container entry point. Its main role
# is to create a user able to execute docker commands from the container
# connecting to the host docker socket at runtime.
#
# The invoker needs to pass the user ID using the variable NXF_USRMAP.
# When this variable is defined, it creates a new user in the container
# with such ID and adds it to the `docker` group, then assigns the docker
# socket file ownership to that user.
#
# Finally it switches the `nextflow` user using the `su` command and
# executes the original target command line.
# 
# authors:
#  Paolo Di Tommaso
#  Emilio Palumbo
#

# enable debugging
[[ "$NXF_DEBUG_ENTRY" ]] && set -x

# wrap cli args with single quote to avoid wildcard expansion
cli=''; for x in "$@"; do cli+="'$x' "; done

# the NXF_USRMAP hold the user ID in the host environment 
if [[ "$NXF_USRMAP" ]]; then
# create a `nextflow` user with the provided ID 
# then change the docker socker ownership to `nextflow` user 
addgroup docker
adduser -u $NXF_USRMAP -G docker -s /bin/bash -D nextflow
chown nextflow /var/run/docker.sock  
# finally run the target command with `nextflow` user
su nextflow << EOF
[[ "$NXF_DEBUG_ENTRY" ]] && set -x
exec bash -c "$cli"
EOF

# otherwise just execute the command
else 
exec bash -c "$cli"
fi