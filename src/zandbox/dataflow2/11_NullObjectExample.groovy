import groovyx.gpars.dataflow.DataflowQueue
import org.codehaus.groovy.runtime.NullObject

/**
 *  --== Null values (skip value) ==--
 *
 *
 *  + If a chained function returns a null value, it is normally passed along the pipeline as a valid value.
 *  + To indicate to the operator that no value should be passed further down the pipeline,
 *    a NullObject.nullObject instance must be returned (in other words it is skipped)
 *
 */

final DataflowQueue queue1 = new DataflowQueue()
final DataflowQueue queue2 = new DataflowQueue()
final odd = {num ->
    if (num == 5) return null  //null values are normally passed on
    if (num % 2 != 0) return num
    else return NullObject.nullObject  //this value gets blocked
}

queue1.chainWith odd into queue2
(1..5).each {queue1 << it}
assert 1 == queue2.val
assert 3 == queue2.val
assert null == queue2.val