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

package nextflow.k8s
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskRun
/**
 * Implements a BASH wrapper for tasks executed by kubernetes cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class K8sWrapperBuilder extends BashWrapperBuilder {

    K8sWrapperBuilder(TaskRun task) {
        super(task)
    }

    /**
     * only for testing purpose -- do not use
     */
    protected K8sWrapperBuilder() {}

    @Override
    boolean fixOwnership() {
        containerConfig.fixOwnership
    }

    @Override
    Map<String,Path> getResolvedInputs() {
        super.getResolvedInputs()
    }


}
