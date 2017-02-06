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