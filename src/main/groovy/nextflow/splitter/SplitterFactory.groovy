/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
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

    @Lazy
    static private Session session = { Global.session as Session } ()

    /**
     * Creates a splitter object by specifying the strategy name
     *
     * @param strategy The splitting strategy e.g. {@code 'fasta'}, {@code 'fastq'}, etc
     * @param object
     *      A map containing named parameters used to initialize the splitter object.
     *      See {@link AbstractSplitter#options(java.util.Map)}
     */
    static SplitterStrategy create( String strategy, String method = null ) {
        assert strategy
        String name = strategy.contains('.') ? strategy : "nextflow.splitter.${strategy.capitalize()}Splitter"
        try {
            def clazz = (Class<SplitterStrategy>) Class.forName(name)
            String operator = method ? method + strategy.capitalize() : null
            create(clazz, operator )
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
        def method = methodName.substring(0,5)
        def qualifier = methodName.substring(5)
        def splitter = create(qualifier, method)
        if( !splitter )
            throw e

        // converts args array to options map
        def opt = argsToOpt(args)

        /*
         * call the 'split'
         */
        if( method == 'split' ) {
            // when the  target obj is a channel use call
            if( obj instanceof DataflowReadChannel ) {
                def outbound = splitOverChannel( obj, splitter, opt )
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
     * Implements dynamic method extension for splitter operators
     *
     * @param source
     * @param splitter
     * @param opts
     * @return
     */
    static protected splitOverChannel( DataflowReadChannel source, SplitterStrategy splitter, Map opts )  {

        def strategy = splitter as AbstractSplitter

        if( opts.containsKey('autoClose') ) throw new IllegalArgumentException('Parameter `autoClose` do not supported')
        if( opts.into && !(opts.into instanceof DataflowQueue) ) throw new IllegalArgumentException('Parameter `into` must reference a channel object')

        // create a new DataflowChannel that will receive the splitter entries
        DataflowQueue resultChannel
        if( opts.into ) {
            resultChannel = opts.into as DataflowQueue
        }
        else {
            resultChannel = opts.into = new DataflowQueue()
        }

        // turn off channel auto-close
        opts.autoClose = false

        // set the splitter strategy options
        strategy.options(opts)

        DataflowHelper.subscribeImpl ( source, [
                onNext: { entry -> strategy.target(entry).apply() },
                onComplete: { resultChannel << Channel.STOP }
        ] )

        // return the resulting channel
        return resultChannel
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
