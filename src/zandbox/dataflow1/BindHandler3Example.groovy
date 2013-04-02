import groovyx.gpars.dataflow.DataflowQueue


/**
 * --== wheneverBound ==--
 *
 * + Dataflow queues and broadcasts also support a wheneverBound method to register
 *   a closure or a message handler to run each time a value is bound to them.
 */

bookingPromise = new DataflowQueue<>()

bookingPromise.wheneverBound {booking -> println "printAgenda: $booking" }
bookingPromise.wheneverBound {booking -> println "sendMeAnEmailTo: $booking" }
bookingPromise.wheneverBound {booking -> println "updateTheCalendar: $booking" }


bookingPromise << 1 << 2

sleep 1000