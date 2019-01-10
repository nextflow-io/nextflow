#!/usr/bin/env nextflow
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

echo true

list1 = [1,2]
list2 = ['Hola', 'Ciao']
list3 = ['alpha','beta','delta']

process hola {

    input:
    val x from list1
    each y from list2
    each z from list3

    """
    echo 'x: $x; y: $y; z: $z'
    """

}

process foo {
    echo true

    input:
    each v from Channel.from([["a","b"],["c","d"]])

    """
    echo foo $v
    """
}
