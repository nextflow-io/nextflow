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