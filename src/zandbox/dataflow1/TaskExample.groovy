/**
 *  --== Dataflow Task ==--
 *
 *  The Dataflow tasks give you an easy-to-grasp abstraction of mutually-independent logical tasks or threads,
 *  which can run concurrently and exchange data solely through Dataflow Variables, Queues, Broadcasts and Streams.
 *
 *  Dataflow tasks with their easy-to-express mutual dependencies and inherently sequential body could also be used as
 *  a practical implementation of UML Activity Diagrams .
 */


import static groovyx.gpars.GParsPool.*
import groovyx.gpars.dataflow.DataflowVariable
import static groovyx.gpars.dataflow.Dataflow.task

/**
 * A simple mashup sample, downloads content of three websites
 * and checks how many of them refer to Groovy.
 */

def dzone = new DataflowVariable()
def jroller = new DataflowVariable()
def theserverside = new DataflowVariable()

task {
    println 'Started downloading from DZone'
    dzone << 'http://www.dzone.com'.toURL().text
    println 'Done downloading from DZone'
}

task {
    println 'Started downloading from JRoller'
    jroller << 'http://www.jroller.com'.toURL().text
    println 'Done downloading from JRoller'
}

task {
    println 'Started downloading from TheServerSide'
    theserverside << 'http://www.theserverside.com'.toURL().text
    println 'Done downloading from TheServerSide'
}

task {
    withPool {
        println "Number of Groovy sites today: " +
                ([dzone, jroller, theserverside].findAllParallel {
                    it.val.toUpperCase().contains 'GROOVY'
                }).size()
    }
}.join()