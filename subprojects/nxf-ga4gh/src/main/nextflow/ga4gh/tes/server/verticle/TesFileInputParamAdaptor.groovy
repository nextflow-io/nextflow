/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.ga4gh.tes.server.verticle

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.script.BaseInParam
import nextflow.script.FileInParam

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TesFileInputParamAdaptor extends FileInParam {

    static private Binding binding0 = new Binding()

    static private List holder0 = []

    String name

    TesFileInputParamAdaptor(String name) {
        super(binding0, holder0)
        this.name = name
    }

    String toString() { name }

    @Override
    String getName() { name }


    @Override
    DataflowReadChannel getInChannel() {
        throw new UnsupportedOperationException()
    }

    @Override
    BaseInParam from(Object obj) {
        throw new UnsupportedOperationException()
    }

    @Override
    protected void lazyInit() {
        throw new UnsupportedOperationException()
    }

    @Override
    BaseInParam bind(Object obj) {
        throw new UnsupportedOperationException()
    }

    @Override
    def decodeInputs(List inputs) {
        throw new UnsupportedOperationException()
    }
}
