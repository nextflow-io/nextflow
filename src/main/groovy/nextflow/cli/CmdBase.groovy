/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cli

import com.beust.jcommander.Parameter

/**
 * Implement command shared methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class CmdBase implements Runnable {

    private Launcher launcher

    abstract String getName()

    Launcher getLauncher() { launcher }

    void setLauncher( Launcher value ) { this.launcher = value }

    @Parameter(names=['-h','-help'], description = 'Print the command usage', arity = 0, help = true)
    boolean help
}