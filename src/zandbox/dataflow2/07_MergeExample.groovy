/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

import javax.xml.crypto.Data

import groovyx.gpars.dataflow.DataflowQueue

/**
 * --== Merging channels ==--
 *
 * + Merging allows you to join multiple read channels as inputs for a single dataflow operator.
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

def maleChannel = new DataflowQueue()
def femaleChannel = new DataflowQueue()
def mortgageCandidatesChannel = new DataflowQueue()

maleChannel.merge(femaleChannel) {m, f -> "${m} and ${f}" }.into(mortgageCandidatesChannel)

maleChannel << 'John'
maleChannel << 'Phillip'

femaleChannel << 'Rose'
femaleChannel << 'Kate'



println "Mortgage: " + mortgageCandidatesChannel.getVal()
println "Mortgage: " + mortgageCandidatesChannel.getVal()