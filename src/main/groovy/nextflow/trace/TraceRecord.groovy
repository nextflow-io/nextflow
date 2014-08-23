/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.trace

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
  * This object represent holds the information of a single process run,
  * its content is saved to a trace file line
  *
  * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
  */
@ToString
@EqualsAndHashCode
class TraceRecord {

    def taskId
    def nativeId
    def status
    def exit
    def hash
    String name
    long submit
    long start
    long complete
    def mem
    def cpu

}