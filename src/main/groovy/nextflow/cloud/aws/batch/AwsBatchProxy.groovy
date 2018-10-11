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

package nextflow.cloud.aws.batch

import com.amazonaws.services.batch.AWSBatch
import nextflow.util.ClientProxyThrottler
import nextflow.util.ThrottlingExecutor
/**
 * Implements a AWS Batch client proxy that handle all API invocations
 * through the provided executor service
 *
 * WARN: the caller class/method should not be compile static
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchProxy extends ClientProxyThrottler<AWSBatch> {

    @Delegate(deprecated=true)
    private AWSBatch target

    AwsBatchProxy(AWSBatch client, ThrottlingExecutor executor) {
        super(client, executor, [describeJobs: 10 as Byte]) // note: use higher priority for `describeJobs` invocations
        this.target = client
    }

}
