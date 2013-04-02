import static groovyx.gpars.dataflow.Dataflow.task

/**
 * --== multiple whenBound ==--
 *
 * + nothing prevents you from having more of such handlers for a single promise:
 *   They will all trigger in parallel once the promise has a concrete value:
 *
 */

bookingPromise = task {
    return [1,2,2,3]
}

bookingPromise.whenBound {booking -> println "printAgenda: $booking" }
bookingPromise.whenBound {booking -> println "sendMeAnEmailTo: $booking" }
bookingPromise.whenBound {booking -> println "updateTheCalendar: $booking" }

bookingPromise >> { println "more1: $it" }
bookingPromise >> { println "done: $it" }


sleep(1000)