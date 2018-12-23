/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
