#!/bin/bash
#
#  Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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