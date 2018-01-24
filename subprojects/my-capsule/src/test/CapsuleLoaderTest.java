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

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CapsuleLoaderTest {

    @Test
    public void testExtendDeps () {
        assertEquals( CapsuleLoader.extendDepsWith(null, null), null);
        assertEquals( CapsuleLoader.extendDepsWith("a:a:1.0", null), Arrays.asList("a:a:1.0"));
        assertEquals( CapsuleLoader.extendDepsWith("a:a:1.0 b:b:2.0", null), Arrays.asList("a:a:1.0","b:b:2.0"));
        assertEquals( CapsuleLoader.extendDepsWith("a:a:1.0   b:b:2.0  ", null), Arrays.asList("a:a:1.0","b:b:2.0"));
        assertEquals( CapsuleLoader.extendDepsWith("x:x:1.0 y:y:2.0", Arrays.asList("a:a:1.0", "b:b:2.0")), Arrays.asList("a:a:1.0","b:b:2.0","x:x:1.0","y:y:2.0"));
    }


    @Test
    public void testExtendClasspath() {
        assertEquals( CapsuleLoader.extendClassPathWith(null, null), null );
        assertEquals( CapsuleLoader.extendClassPathWith("x", null), Arrays.asList("x"));
        assertEquals( CapsuleLoader.extendClassPathWith("x:y ", null), Arrays.asList("x","y"));
        assertEquals( CapsuleLoader.extendClassPathWith("x :: y ", null), Arrays.asList("x","y"));
        assertEquals( CapsuleLoader.extendClassPathWith("x:y:z", Arrays.asList("lib1.jar","lib2.jar")), Arrays.asList( "lib1.jar","lib2.jar","x","y","z"));

    }

}
