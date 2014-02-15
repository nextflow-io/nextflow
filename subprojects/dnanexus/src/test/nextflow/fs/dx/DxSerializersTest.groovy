package nextflow.fs.dx

import nextflow.fs.dx.api.DxApi
import nextflow.util.KryoHelper
import spock.lang.Specification
/*
 * Copyright (c) 2012, the authors.
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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DxSerializersTest extends Specification {


    def testSerialization() {

        setup:
        KryoHelper.register(DxPath, DxPathSerializer)
        KryoHelper.register(DxFileSystem, DxFileSystemSerializer)

        //provider.fileSystems =
        DxFileSystemProvider provider = new DxFileSystemProvider('123', Mock(DxApi))
        DxFileSystemProvider.instance = provider

        def fs = new DxFileSystem(provider,'123','sandbox')
        def pathObj = new DxPath(fs, '/home/hola')

        def store = File.createTempFile('dxser',null)
        store.deleteOnExit()

        when:
        KryoHelper.write(pathObj, store)
        def copy = (DxPath)KryoHelper.read(store)

        then:
        copy == pathObj
        pathObj.getFileName().toString() == 'hola'
        pathObj.toString() == 'dxfs://sandbox:/home/hola'


        cleanup:
        DxFileSystemProvider.instance = null


    }


}
