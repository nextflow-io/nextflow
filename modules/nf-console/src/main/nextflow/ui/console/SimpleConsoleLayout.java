/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
