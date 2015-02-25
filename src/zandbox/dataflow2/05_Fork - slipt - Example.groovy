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

import static groovyx.gpars.dataflow.Dataflow.task

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel

/**
 *  --== Forking the data flow ==--
 *
 *  + When a need comes to copy the output of a pipeline/channel into more than one following pipeline/channel,
 *    the split() method will help you:
 *
 */

def toUpperCase = {s -> s.toUpperCase()}

final encrypt = new DataflowQueue()
final DataflowWriteChannel messagesToSave = new DataflowQueue()
final DataflowWriteChannel messagesToLog = new DataflowQueue()
encrypt.chainWith toUpperCase split(messagesToSave, messagesToLog)

encrypt << "I need to keep this message secret!"
encrypt << "GPars can build operator pipelines really easy"


println 'log: ' + messagesToLog.getVal() + ' - '+ messagesToLog.getVal()
println 'save: ' + messagesToSave.getVal() + ' - '+ messagesToSave.getVal()


