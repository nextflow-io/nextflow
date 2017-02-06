/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import groovyx.gpars.dataflow.DataflowVariable


/**
 *  --== Bind handler ==--
 *
 *  + Bind handlers can be registered on all dataflow channels (variables, queues or broadcasts) either
 *    using the >> operator and the then() or the whenBound() methods.
 *
 *  + They will be run once a value is bound to the variable.
 */

def a = new DataflowVariable()
a >> {println "The variable has just been bound to $it"}
a.whenBound {println "Just to confirm that the variable has been really set to $it"}

a << 1

sleep 1000
