/*
 * Copyright 2013-2023, Seqera Labs
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

package io.seqera.wave.config

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SpackOptsTest extends Specification {

    def 'check spack default options' () {
        given:
        def opts = new SpackOpts()
        expect:
        opts.checksum
        opts.builderImage == SpackOpts.DEFAULT_SPACK_BUILDER_IMAGE
        opts.runnerImage == SpackOpts.DEFAULT_SPACK_RUNNER_IMAGE
        opts.osPackages == SpackOpts.DEFAULT_SPACK_OSPACKAGES
        opts.commands == null
    }

    def 'check spack custom opts' () {
        given:
        def opts = new SpackOpts([
                checksum:false,
                builderImage: 'my/builder:image',
                runnerImage: 'my/runner:image',
                osPackages:  'my-os-packages',
                commands: ['run','--this','--that']
        ])

        expect:
        !opts.checksum
        and:
        opts.builderImage == 'my/builder:image'
        opts.runnerImage == 'my/runner:image'
        opts.osPackages == 'my-os-packages'
        and:
        opts.commands == ['run','--this','--that']
    }
}
