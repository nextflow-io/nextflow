/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.k8s.model

import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a K8s Secret file mount
 *
 * https://kubernetes.io/docs/concepts/configuration/secret/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class PodMountSecret {

    String mountPath

    String fileName

    String secretName

    String secretKey

    PodMountSecret(String secret, String mount) {
        assert secret
        assert mount

        final path = Paths.get(mount)
        final tokens = secret.tokenize('/')
        secretName = tokens[0].trim()
        secretKey = tokens.size()>1 ? tokens[1].trim() : null
        if( secretKey ) {
            mountPath = path.parent.toString()
            fileName = path.fileName.toString()
        }
        else {
            mountPath = path.toString()
        }
    }

    PodMountSecret(Map entry) {
        this(entry.secret as String, entry.mountPath as String)
    }

}
