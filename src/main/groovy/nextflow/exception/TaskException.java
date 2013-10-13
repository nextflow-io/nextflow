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

package nextflow.exception;

/**
 * Reports a generic error during the 'Task' validation step
 *
 * Note: THIS IS A PLAIN JAVA CLASS due to this bug
 * http://blog.proxerd.pl/article/how-to-fix-incompatibleclasschangeerror-for-your-groovy-projects-running-on-jdk7
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class TaskException extends Exception {

    public TaskException(String message) {
        super(message);
    }

    public TaskException(String message, Throwable cause) {
        super(message,cause);
    }
}

