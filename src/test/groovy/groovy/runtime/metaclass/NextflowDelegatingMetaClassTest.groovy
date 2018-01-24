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

package groovy.runtime.metaclass
import java.nio.file.Files

import nextflow.Channel
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowDelegatingMetaClassTest extends Specification {


    /*
     * 'isEmpty' method is implemented by the NextflowDelegatingMetaClass class
     */
    def testEmptyAndIsEmpty() {

        when:
        def file = File.createTempFile('hello','file')
        then:
        file.empty()
        file.isEmpty()
        when:
        file.text = 'hello'
        then:
        !file.empty()
        !file.isEmpty()

        when:
        def path = Files.createTempFile('hello','path')
        then:
        path.empty()
        path.isEmpty()
        when:
        path.text = 'hello'
        then:
        !path.empty()
        !path.isEmpty()

        cleanup:
        file?.delete()
        path?.delete()
    }


    def testStringSplitText() {

        when:
        def str = "hello\nworld"
        def list = str.splitText(into: []) { it.trim() }
        then:
        list == ['hello','world']

    }

    def testStringSplitEach() {

        when:
        def str = "hello\nworld"
        def list = str.splitText(into: [], each: { it.trim().reverse() })
        then:
        list == ['olleh','dlrow']

    }

    def testGStringSplitText() {

        when:
        def x = 'hello'
        def y = 'world'
        def str = "$x\n$y"
        def list = str.splitText(into: []) { it.trim() }
        then:
        list == ['hello','world']

    }


    def testFileSplitText() {

        given:
        def file = File.createTempFile('testSplit', null)

        when:
        file.text = 'Hello\nworld\n!'
        def list = file.splitText(into: []) { it.trim() }
        then:
        list == ['Hello','world','!']

        when:
        file.text = 'Hello\nworld\n!'
        list = file.splitText(into: [], by:2)
        then:
        list == ['Hello\nworld\n','!\n']

        when:
        file.text = 'Hello\nworld\n!'
        list = file.splitText(into: [], each: { it.trim().reverse() })
        then:
        list == ['olleH','dlrow','!']

        cleanup:
        file?.delete()

    }


    def testPathSplitText() {

        given:
        def path = Files.createTempFile('splitPath', null)

        when:
        path.text = 'Hello\nworld\n!'
        def list = path.splitText(into: [])
        then:
        list == ['Hello\n','world\n','!\n']

        when:
        path.text = 'Hello\nworld\n!'
        def list2 = path.splitText(into: [], by:2 )
        then:
        list2 == ['Hello\nworld\n','!\n']

        when:
        path.text = 'Hello\nworld\n!'
        def list3 = path.splitText(into: [], each: { it.trim().reverse() })
        then:
        list3 == ['olleH','dlrow','!']

        cleanup:
        path?.delete()

    }

    def testCountText() {

        when:
        def str = '''
            line 1
            line 2
            line 3
            line 4
            line 5
            '''
            .stripIndent().strip()

        then:
        str.countText() == 5
        str.countText(by:2) == 3


        when:
        def str2 = '''
            line 6
            line 7
            line 8
            '''
                .stripIndent().strip()

        def str3 = '''
            line 9
            line 10
            line 11
            '''
                .stripIndent().strip()

        def result = Channel.from(str, str2, str3).countText()
        then:
        result.val == 11


    }

    def testCountFasta() {

        when:
        def str = '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            >1ihvA
            NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
            '''
            .stripIndent()

        then:
        str.countFasta() == 5
        str.countFasta(by:4) == 2



        when:
        def str2 = '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            '''
            .stripIndent()

        def result = Channel.from( str, str2 ).countFasta()
        then:
        result.val == 8


    }


}
