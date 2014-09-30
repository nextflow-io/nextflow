package test
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
import org.iq80.leveldb.*;
import static org.iq80.leveldb.impl.Iq80DBFactory.*;
import java.io.*;

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

Options options = new Options();
options.createIfMissing(true);
DB db = factory.open(new File("example"), options);
try {
    db.put(bytes("alpha"), bytes("rocks"));
    db.put(bytes("beta"), bytes("pops"));
    db.put(bytes("delta"), bytes("metal"));

    DBIterator iterator = db.iterator();
    try {
        for(iterator.seekToFirst(); iterator.hasNext(); iterator.next()) {
            String key = asString(iterator.peekNext().getKey());
            String value = asString(iterator.peekNext().getValue());
            System.out.println(key+" = "+value);
        }
    } finally {
        // Make sure you close the iterator to avoid resource leaks.
        iterator.close();
    }

}
finally {
    // Make sure you close the db to shutdown the
    // database and avoid resource leaks.
    db.close();
}