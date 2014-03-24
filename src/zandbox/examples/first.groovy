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

import nextflow.Channel

Channel.from( 1, 2, 3 ).first() .subscribe { println it }

//
Channel.from( 1, 2, 'a', 'b', 3 ) .first( String ) .subscribe { println it }

//
Channel.from( 'a', 'aaa', 'aaa' ) .first( ~/aa.*/) .subscribe { println it }

Channel.from( 1,2,3,4,5 ) .first { it > 3 } .subscribe { println it }


sleep 100