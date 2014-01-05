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

package nextflow.script

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileSpec {

    private String fFilePattern

    private String fType = 'file'

    private boolean fCreate

    FileSpec filePattern( String pattern ) {
        this.fFilePattern = pattern
        return this
    }

    String getFilePattern() {
        return fFilePattern
    }

    boolean getCreate() {
        return fCreate
    }

    boolean isDirectory() {
        return fType == 'dir'
    }

    boolean isFile() {
        return fType == 'file'
    }

    FileSpec create(boolean flag) {
        this.fCreate = flag
        return this
    }

    FileSpec type(String value) {
        assert value in ['file','dir']
        fType = value
        return this
    }


}
