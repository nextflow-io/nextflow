/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.trace

import java.nio.file.Path

import groovy.json.JsonBuilder
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Render the preview report when running a pipeline
 * in preview mode.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PreviewReportWriter implements TraceObserver {

    static private final List<String> DEF_DIRECTIVES = ['container', 'cpus', 'memory', 'time']

    static private final String DEF_FILE_NAME = "preview-${TraceHelper.launchTimestampFmt()}.json"

    private DAG dag

    private List<String> directives

    private Path file

    private boolean overwrite

    static PreviewReportWriter create(Session session) {
        final directives = session.config.navigate('preview.directives', DEF_DIRECTIVES) as List
        def file = session.config.navigate('preview.file', DEF_FILE_NAME)
        file = (file as Path).complete()
        final overwrite = session.config.navigate('preview.overwrite', false) as boolean

        return new PreviewReportWriter(session.dag, directives, file, overwrite)
    }

    PreviewReportWriter( DAG dag, List<String> directives, Path file, boolean overwrite ) {
        this.dag = dag
        this.directives = directives
        this.file = file
        this.overwrite = overwrite

        // check file existance
        final attrs = FileHelper.readAttributes(file)
        if( attrs ) {
            if( overwrite && (attrs.isDirectory() || !file.delete()) )
                throw new AbortOperationException("Unable to overwrite existing preview file: ${file.toUriString()}")
            else if( !overwrite )
                throw new AbortOperationException("Preview file already exists: ${file.toUriString()} -- enable `preview.overwrite` in your config file to overwrite existing preview files")
        }
    }

    void render() {
        // get preview data
        final jsonData = [:]

        for( def vertex : dag.vertices ) {
            // skip nodes that are not processes
            final process = vertex.process
            if( !process )
                continue

            // get preview task config
            final taskConfig = process.getPreviewConfig()

            // get preview for each directive
            def config = [:]

            for( def directive : directives ) {
                // try to resolve directive value
                def value = null
                try {
                    value = taskConfig.get(directive)
                }
                catch( Exception e ) {
                    log.warn1 "Unable to preview directive `${directive}` for process `${process.name}`: ${e}"
                }

                // convert custom types to string values
                if( value instanceof Duration || value instanceof MemoryUnit )
                    value = value.toString()

                // add directive value if it was resolved
                if( value != null )
                    config[directive] = value
            }

            jsonData[process.name] = config
        }

        // write preview data to file
        file.text = new JsonBuilder(jsonData).toString()
    }

}
