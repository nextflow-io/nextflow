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

package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.extension.DataflowExtensions

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SplitterFactory {

    private static final String SPLIT_PREFIX = 'split'

    /**
     *
     * @param strategy The splitting strategy e.g. {@code 'fasta'}, {@code 'fastq'}, etc
     * @param object
     */
    @Memoized
    static SplitterStrategy create( String strategy, Map options = [:] ) {
        assert strategy
        String name = strategy.contains('.') ? strategy : "nextflow.splitter.${strategy.capitalize()}Splitter"
        try {
            def clazz = (Class<SplitterStrategy>) Class.forName(name)
            create(clazz, options)
        }
        catch( ClassNotFoundException e ) {
            log.debug "Cannot find any class implementing a split strategy for: $strategy"
            return null
        }

    }

    @Memoized
    static SplitterStrategy create( Class<? extends SplitterStrategy> strategy, Map options = [:] ) {
        (SplitterStrategy) strategy.newInstance( [options] as Object[] )
    }

    static tryFallbackWithSplitter( obj, String methodName, Object[] args, Exception e ) {

        // verifies that is a splitter method and get splitter qualifier
        if( !methodName.startsWith(SPLIT_PREFIX) || !(methodName.size()>SPLIT_PREFIX.length()) )
            throw e

        // load the splitter class
        def qualifier = methodName.substring(SPLIT_PREFIX.length())
        def splitter = create(qualifier)
        if( !splitter )
            throw e

        // converts args array to options map
        def opt = argsToOpt(args)

        // when the  target obj is a channel use call
        if( obj instanceof DataflowReadChannel ) {
            splitOverChannel( obj, splitter, opt )
        }
        else {
            // invokes the splitter
            splitter.options(opt) .target(obj) .split()
        }

    }

    static protected splitOverChannel( DataflowReadChannel source, SplitterStrategy splitter, Map opt )  {

        def strategy = splitter as AbstractSplitter

        // create a new DataflowChannel that will receive the splitter entries
        DataflowQueue resultChannel = opt.into = new DataflowQueue()

        // turn off channel auto-close
        opt.autoClose = false

        // set the splitter strategy options
        strategy.options(opt)

        int count = 0
        def splitEntry = { entry ->
            def obj = strategy.normalizeType(entry)
            strategy.apply(obj, count)
        }

        DataflowExtensions.subscribe ( source, [onNext: splitEntry, onComplete: { resultChannel << Channel.STOP }] )

        // return the resulting channel
        return resultChannel
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
