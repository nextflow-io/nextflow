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

package nextflow.executor

import groovy.transform.EqualsAndHashCode
import nextflow.Session

/**
 * Wrap the client session and classpath data, in order to enable
 * remote closure execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
class GgClient implements Serializable {

    private static final long serialVersionUID = - 3956143328835315200L ;

    final UUID sessionId

    final List<URL> classpath

    GgClient() { }

    GgClient(Session session) {
        sessionId = session.getUniqueId()
        classpath = session.getClasspath()
    }

    /** only for test */
    protected GgClient( UUID uuid, List<URL> classpath  ) {
        this.sessionId = uuid
        this.classpath = classpath
    }

}