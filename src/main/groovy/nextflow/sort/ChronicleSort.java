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

package nextflow.sort;

import java.io.File;
import java.io.IOException;

import net.openhft.chronicle.map.ChronicleMap;
import net.openhft.chronicle.map.ChronicleMapBuilder;
import nextflow.extension.FileEx;
import nextflow.util.KryoHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implements external sort sorting values into a Chronicle persistent map
 * https://github.com/OpenHFT/Chronicle-Map
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChronicleSort<V> extends BigSort<V> {

    private static final Logger log = LoggerFactory.getLogger(ChronicleSort.class);

    private ChronicleMap<Long,byte[]> data;

    private long entries;

    /**
     * Define the number of expected entries to be stored in the map
     *
     * @param value The number of expected entries
     * @return The sort object itself
     */
    public ChronicleSort<V> entries( long value ) {
        entries = value;
        return this;
    }

    /**
     * Create the underlying Chronicle map persisted structure
     * @return
     * @throws IOException
     */
    public ChronicleSort<V> create() throws IOException {
        super.create();

        /*
         * creates the chronicle-map storage
         */
        new File(getTempDir(), "data").mkdirs();
        File file = new File(getTempDir(), "data/chronicle.dat");
        ChronicleMapBuilder<Long,byte[]> builder = ChronicleMapBuilder.of(Long.class, byte[].class);
        if( entries > 0 ) {
            log.debug("Chronicle entries: {}", entries);
            builder.entries(entries);
        }

        data = builder.create(file);

        return this;
    }

    @Override
    protected int put(long key, V value) {
        byte[] bytes = KryoHelper.serialize(value);
        data.put(key,bytes);
        return 8 + bytes.length;
    }

    @Override
    @SuppressWarnings("unchecked")
    protected V get(long key) {
        byte[] raw = data.get(key);
        return (V)KryoHelper.deserialize(raw);
    }

    /**
     * Close underlying structures
     */
    public void close() {
        log.debug("ChronicleSort#close");
        FileEx.closeQuietly(data);
        log.debug("ChronicleSort#close -- map closed");
        super.close();
        log.debug("ChronicleSort#close -- file deleted");
    }
}