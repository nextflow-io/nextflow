/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package com.upplication.s3fs.experiment;

import java.math.BigInteger;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import com.google.common.annotations.Beta;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Tentative for transfer rate meter
 *
 * DON'T USE, IT ADDS NOT NEGLIGIBLE EXECUTION OVERHEAD
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Beta
public class TransferRateMeter {

    static private final Logger log = LoggerFactory.getLogger(TransferRateMeter.class);

    static private AtomicBigInteger accumulator = new AtomicBigInteger();
    static long UPDATE_DELAY_MS = 5_000;
    static long timeLastPrint = System.currentTimeMillis();
    static BigInteger valueLastPrint = BigInteger.ZERO;

    static final private DecimalFormatSymbols formatSymbols;
    final static public String[] UNITS = { "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB" };

    static {
        formatSymbols = new DecimalFormatSymbols();
        formatSymbols.setDecimalSeparator('.');
    }


    private long count;
    private long lastMillis;

    public void inc(long value) {
        final long now = System.currentTimeMillis();
        count += value;
        if( now-lastMillis > UPDATE_DELAY_MS) {
            accumulate(count);
            count=0;
            lastMillis = now;
        }
    }

    public void inc() {
        inc(1);
    }

    static private BigInteger accumulate( long value ) {
        accumulator.getAndIncrement(value);
        long now = System.currentTimeMillis();
        if( now- timeLastPrint < UPDATE_DELAY_MS )
            return null;
        synchronized (TransferRateMeter.class) {
            long ts = System.currentTimeMillis();
            long timeDelta = ts - timeLastPrint;
            BigInteger current = accumulator.get();
            long divisor = timeDelta / 1000;
            if( divisor > 0 ) {
                BigInteger rate = current.subtract(valueLastPrint).divide( BigInteger.valueOf(divisor) );
                log.debug("Download transfer rate: {}", toUnitString(rate.longValue()));
                timeLastPrint = ts;
                valueLastPrint = current;
                return rate;
            }
            return null;
        }
    }

    static private String toUnitString(long size) {
        if(size <= 0) {
            return "0";
        }
        // see http://stackoverflow.com/questions/2510434/format-bytes-to-kilobytes-megabytes-gigabytes
        int digitGroups = (int) (Math.log10(size) / Math.log10(1024));
        DecimalFormat formatter = new DecimalFormat("0.#", formatSymbols);
        formatter.setGroupingUsed(false);
        return formatter.format(size / Math.pow(1024, digitGroups)) + " " + UNITS[digitGroups];
    }
}
