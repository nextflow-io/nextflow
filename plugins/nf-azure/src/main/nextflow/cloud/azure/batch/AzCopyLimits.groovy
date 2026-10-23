/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.cloud.azure.batch

import java.math.RoundingMode

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.MemoryUnit

/**
 * Sizes azcopy for the resources of a single task.
 *
 * By default azcopy sizes its buffer from the CPUs of the host, ignoring the
 * memory limit of the task container it runs in, and gets OOM-killed while
 * staging files for tasks with little memory.
 *
 * @author Jonathan Manning <jonathan.manning@seqera.io>
 */
@CompileStatic
@EqualsAndHashCode
@ToString(includeNames = true, includePackage = false)
class AzCopyLimits {

    // the page cache of the files being written is charged to the container too
    static final long RESERVED_MB = 128
    static final long PROCESS_OVERHEAD_MB = 64
    static final long MIN_BUFFER_MB = 32
    // azcopy waits indefinitely for buffer space when a block does not fit
    static final int MIN_BLOCKS_IN_BUFFER = 4
    static final long MIN_TRANSFER_MB = 192

    final int parallelTransfers
    final long bufferMb

    AzCopyLimits(int parallelTransfers, long bufferMb) {
        this.parallelTransfers = parallelTransfers
        this.bufferMb = bufferMb
    }

    /**
     * @param memory The task memory, or {@code null} when not specified
     * @param maxParallelTransfers The configured maximum number of parallel transfers
     * @param blockSize The azcopy block size in MB
     * @return The azcopy limits for the task, or {@code null} when the task memory is not specified
     */
    static AzCopyLimits of(MemoryUnit memory, int maxParallelTransfers, String blockSize) {
        if( !memory )
            return null
        final blockSizeMb = blockSize?.isNumber() ? (long)Math.ceil(blockSize as double) : 0L
        final budget = memory.toMega().intdiv(2) - RESERVED_MB
        final parallel = budget > 0
                ? Math.max(1, Math.min(maxParallelTransfers, (int)(budget.intdiv(MIN_TRANSFER_MB))))
                : 1
        final minBuffer = Math.max(MIN_BUFFER_MB, MIN_BLOCKS_IN_BUFFER * blockSizeMb)
        final buffer = budget > 0
                ? Math.max(minBuffer, budget.intdiv(parallel) - PROCESS_OVERHEAD_MB)
                : minBuffer
        return new AzCopyLimits(parallel, buffer)
    }

    /**
     * @return The buffer size in the format expected by the {@code AZCOPY_BUFFER_GB} variable
     */
    String getBufferGb() {
        return new BigDecimal(bufferMb).divide(new BigDecimal(1024), 3, RoundingMode.UP).stripTrailingZeros().toPlainString()
    }
}
