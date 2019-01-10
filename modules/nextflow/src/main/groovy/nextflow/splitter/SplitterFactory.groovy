/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.extension.DataflowHelper
/**
 * Factory class for splitter objects
 *
 * {@link groovy.runtime.metaclass.NextflowDelegatingMetaClass}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SplitterFactory {

    /**
     * Creates a splitter object by specifying the strategy name
     *
     * @param strategy The splitting strategy e.g. {@code 'fasta'}, {@code 'fastq'}, etc
     * @param method The actual method name used to invoke the splitter function
     */
    static SplitterStrategy create( String strategy ) {
        assert strategy

        String name = strategy.startsWith('split') || strategy.startsWith('count') ? strategy = strategy.substring(5) : strategy

        if( !name.contains('.') )
            name = "nextflow.splitter.${strategy.capitalize()}Splitter"

        try {
            def clazz = (Class<SplitterStrategy>) Class.forName(name)
            create(clazz, strategy.capitalize())
        }
        catch( ClassNotFoundException e ) {
            log.debug "Cannot find any class implementing a split strategy for: $strategy"
            return null
        }

    }

    /**
     * Creates a splitter object by specifying the splitter class
     *
     * @param strategy A class implementing {@link SplitterStrategy}
     * @param object
     *      A map containing named parameters used to initialize the splitter object.
     *      See {@link AbstractSplitter#options(java.util.Map)}
     * @return The splitter instance
     */
    static SplitterStrategy create( Class<? extends SplitterStrategy> strategy, String name ) {
        (SplitterStrategy) ( name
                ? strategy.newInstance(name)
                : strategy.newInstance())
    }


    /**
     *  Implements dynamic method extension for counter operators
     *
     * @param source
     * @param splitter
     * @param opt
     * @return
     */
    @PackageScope
    static countOverChannel( DataflowReadChannel source, SplitterStrategy splitter, Map opt )  {

        // create a new DataflowChannel that will receive the splitter entries
        DataflowVariable result = new DataflowVariable ()

        def strategy = splitter as AbstractSplitter

        // set the splitter strategy options
        long count = 0
        if( opt == null ) opt = [:]
        opt.each = { count++ }
        strategy.options(opt)

        def events = new HashMap(2)
        events.onNext = { entry -> strategy.target(entry).apply() }
        events.onComplete = { result.bind(count) }

        DataflowHelper.subscribeImpl ( source, events )

        // return the resulting channel
        return result
    }


}
