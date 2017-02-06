/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package test

import java.nio.file.Paths

import nextflow.sort.LevelDbSort

def dbFolder = Paths.get("/Users/pditommaso/Downloads/db/")
dbFolder.deleteDir()

final collector = new LevelDbSort<String>()
            .tempDir(dbFolder)
            .deleteTempFilesOnClose(false)
            .create()

println("Begin");
final Long now = System.currentTimeMillis();
File text = new File("/Users/pditommaso/Downloads/hs_ref_GRCh38_chr1.fa");// hs_ref_GRCh38_chr1.fa
text.eachLine { String it ->
    collector.add(it)
}
println("End adding - time: ${(System.currentTimeMillis() - now) / 1000} secs -- Stats: ${collector.stats()}" );

final PrintWriter sortFile = new PrintWriter(new FileWriter(dbFolder.resolve("sort.txt").toFile()));
int count=0;
collector.sort {
    count++
    sortFile.println(it)
}
sortFile.close();
println("End sort - time: ${(System.currentTimeMillis() - now) / 1000} secs -- processed: $count items" );

collector.close();

