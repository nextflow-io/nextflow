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

package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.dag.NodeMarker
import nextflow.extension.DataflowHelper
import nextflow.extension.SplitOp
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
     * This method try to invoke a splitter or a counter method dynamically
     *
     * @param obj The target object
     * @param methodName
     *      The splitter or counter method to be invoked. It must start with splitXxx or countXxx
     *      where Xxx represents the splitting format/strategy. For example {@code splitFasta} will invoke
     *      the split method by using the {@link FastaSplitter} class
     * @param args Splitter arguments. See {@link SplitterStrategy#options(java.util.Map)}
     * @param e Exception object to raise if the method is not available
     * @return the splitter result
     */
    static tryFallbackWithSplitter( obj, String methodName, Object[] args, Exception e ) {

        // verifies that is a splitter method and get splitter qualifier
        if( !methodName.startsWith('split') && !methodName.startsWith('count') )
            throw e

        if( methodName.size() == 5 )
            throw e

        // load the splitter class
        def splitter = create(methodName)
        if( !splitter )
            throw e

        // converts args array to options map
        def opt = argsToOpt(args)

        /*
         * call the 'split'
         */
        if( methodName.startsWith('split') ) {
            // when the  target obj is a channel use call
            if( obj instanceof DataflowReadChannel ) {
                def outbound = new SplitOp( obj, methodName, opt ).apply()
                NodeMarker.addOperatorNode(methodName, obj, outbound)
                return outbound
            }
            // invokes the splitter
            else {
                splitter.options(opt) .target(obj) .split()
            }
        }

        /*
         * otherwise handle 'count'
         */
        else {
            // DEPRECATED TO BE REMOVED
            log.warn "Method `$methodName` has been deprecated and it will be removed in a future release"
            // when the  target obj is a channel use call
            if( obj instanceof DataflowReadChannel ) {
                def outbound = countOverChannel( obj, splitter, opt )
                NodeMarker.addOperatorNode(methodName, obj, outbound)
                return outbound
            }
            // invokes the splitter
            else {
                splitter.options(opt) .target(obj) .count()
            }

        }

    }



    /**
     *  Implements dynamic method extension for counter operators
     *
     * @param source
     * @param splitter
     * @param opt
     * @return
     */
    static protected countOverChannel( DataflowReadChannel source, SplitterStrategy splitter, Map opt )  {


        // create a new DataflowChannel that will receive the splitter entries
        DataflowVariable result = new DataflowVariable ()

        def strategy = splitter as AbstractSplitter

        // set the splitter strategy options
        long count = 0
        if( opt == null ) opt = [:]
        opt.each = { count++ }
        strategy.options(opt)

        def splitEntry = { entry ->
            strategy.target(entry).apply()
        }

        DataflowHelper.subscribeImpl ( source, [onNext: splitEntry, onComplete: { result.bind(count) }] )

        // return the resulting channel
        return result
    }



    static protected Map argsToOpt( Object[] args ) {

        Closure closure = null
        Map opt = null

        if( args.size() == 1 ) {
            if( args[0] instanceof Closure )
                closure = args[0] as Closure

            else if( args[0] instanceof Map )
                opt = args[0] as Map

            else
                throw new IllegalArgumentException()
        }
        else if( args.size() == 2 ) {
            opt = args[0] as Map
            closure = args[1] as Closure
        }
        else if( args.size()>2 )
            throw new IllegalArgumentException()

        if( opt == null )
            opt = [:]

        if( closure )
            opt.each = closure

        return opt
    }


}
