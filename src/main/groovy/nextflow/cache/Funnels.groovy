/*
 * Copyright (c) 2012, the authors.
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

package nextflow.cache

import com.google.common.hash.Funnel
import com.google.common.hash.PrimitiveSink

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MapFunnel implements Funnel<Map> {

    @Override
    void funnel(Map map, PrimitiveSink into) {

        map.each { name, value ->
            switch( value?.getClass() ) {
                case CharSequence:
                    into.putString(value as CharSequence)
                    break

                case Integer:
                    into.putInt(value as Integer)
                    break

                case Long:
                    into.putLong(value as Long)
                    break

                case Boolean:
                    into.putBoolean(value as Boolean)
                    break

                case File:
                    into.put

            }
        }

    }
}

class FileFunnel implements Funnel<File> {

    @Override
    void funnel(File file, PrimitiveSink into) {
        assert file

        into
            .putString( file.absolutePath )
            .putLong( file.lastModified() )
            .putLong( file.length() )
    }

}
