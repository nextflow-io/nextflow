/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.processor.tip

import groovy.transform.CompileStatic
import nextflow.processor.TaskRun

/**
 * Implements the default tip provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultTaskTipProvider implements TaskTipProvider {

    static final private List<String> TIPS = [
        'when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line',
        "you can try to figure out what's wrong by changing to the process work dir and showing the script file named `${TaskRun.CMD_SCRIPT}`" as String,
        "view the complete command output by changing to the process work dir and entering the command `cat ${TaskRun.CMD_OUTFILE}`" as String,
        "you can replicate the issue by changing to the process work dir and entering the command `bash ${TaskRun.CMD_RUN}`" as String
    ]

    static final private Random rnd = new Random()

    boolean enabled() { return true }

    @Override
    String suggestTip(List<String> context) {
        return  TIPS[ rnd.nextInt( TIPS.size() ) ]
    }

}
