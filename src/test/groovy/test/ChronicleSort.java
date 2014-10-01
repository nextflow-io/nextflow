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

package test;

import java.io.File;
import java.io.IOException;

import net.openhft.chronicle.map.ChronicleMap;
import net.openhft.chronicle.map.ChronicleMapBuilder;
import nextflow.extension.FileEx;
import nextflow.sort.BigSort;
import nextflow.util.KryoHelper;

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ChronicleSort<V> extends BigSort<V> {

    private ChronicleMap<Long,byte[]> data;

    public BigSort create() throws IOException {
        super.create();
        // creates the chronicle-map storage
        new File(getTempDir(), "data").mkdirs();
        File file = new File(getTempDir(), "data/chronicle.dat");
        data = ChronicleMapBuilder.of(Long.class, byte[].class).file(file).entries(4_000_000).create();

        return this;
    }

    @Override
    protected int put(long key, V value) {
        byte[] bytes = KryoHelper.serialize(value);
        data.put(key,bytes);
        return 8 + bytes.length;
    }

    @Override
    protected V get(long key) {
        byte[] raw = data.get(key);
        return (V)KryoHelper.deserialize(raw);
    }

    public void close() {
        FileEx.closeQuietly(data);
        super.close();
    }
}
