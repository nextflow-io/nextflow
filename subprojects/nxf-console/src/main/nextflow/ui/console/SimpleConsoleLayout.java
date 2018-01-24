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

package nextflow.ui.console;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.CoreConstants;
import ch.qos.logback.core.LayoutBase;

/**
 * Simplified logger layout used to output in the console window area
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class SimpleConsoleLayout extends LayoutBase<ILoggingEvent> {

    @Override
    public String doLayout(ILoggingEvent event) {

        StringBuilder buffer = new StringBuilder(128);
        if( event.getLevel() == Level.INFO ) {
            buffer .append(event.getFormattedMessage()) .append(CoreConstants.LINE_SEPARATOR);
        }

        else {
            buffer
                .append( event.getLevel().toString() ) .append(": ")
                .append(event.getFormattedMessage())
                .append(CoreConstants.LINE_SEPARATOR);
        }

        return buffer.toString();

    }
}
