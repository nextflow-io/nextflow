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

package nextflow.executor
import java.nio.file.Path

import groovy.transform.PackageScope
import nextflow.processor.TaskRun
/**
 * An executor specialised for BSC cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BscExecutor extends LsfExecutor {

    private static String SPECIAL_CHARS = ' []|&!<>'

    /**
     * Defines the job directives that must pre-prepend the job submit script e.g.
     * <pre>
     *  #BSUB -cwd /scratch
     *  #BSUB -q bsc_ls
     *  #BSUB -n 2
     *  #BSUB -R span[hosts=1]
     *  #BSUB -W 01:30
     *  #BSUB -M 4096
     * </pre>
     *
     * @param task The {@link TaskRun} instance for which return the headers script
     * @return A string of new-line separated directives
     */
    String getHeaders(TaskRun task) {

        def result = new StringBuilder()
        def directives = getDirectives(task)

        def len = directives.size()-1
        for( int i=0; i<len; i+=2) {
            def key = directives[i]
            def val = directives[i+1].toString()
            if( val.indexOf(' ') ||  val.indexOf('!') || val.indexOf('[') || val.indexOf(']') )
            result << '#BSUB ' << key << ' ' << val << '\n'
        }

        result.toString()
    }

    /**
     * The command line to submit this job
     * @param task The {@link TaskRun} instance to submit for execution to the cluster
     * @param scriptFile The file containing the job launcher script
     * @return A list representing the submit command line
     */
    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile ) { ['bsub'] }

    /**
     * @return {@code true} since BSC grid requires the script to be piped to the {@code bsub} command
     */
    @Override
    boolean pipeLauncherScript() { true }

    @PackageScope
    String wrap( String str ) {
        boolean needWrap = false

        if( !str ) return str

        for( int i=0; i<SPECIAL_CHARS.size(); i++ ) {
            for( int x=0; x<str.size(); x++ ) {
                if( str.charAt(x) == SPECIAL_CHARS.charAt(i) ) {
                    needWrap = true
                    break
                }
            }
        }

        return needWrap ? "\"$str\"" : str
    }

}
