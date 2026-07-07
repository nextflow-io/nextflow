/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.nix

import spock.lang.Specification
import spock.lang.Unroll

class NixConfigTest extends Specification {

    @Unroll
    def 'should check enabled flag'() {
        given:
        def nix = new NixConfig(CONFIG, ENV)
        expect:
        nix.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]
        and:
        false       | [:]               | [NXF_NIX_ENABLED: 'false']
        true        | [:]               | [NXF_NIX_ENABLED: 'true']
        false       | [enabled: false]  | [NXF_NIX_ENABLED: 'true']  // <-- config has priority
        true        | [enabled: true]   | [NXF_NIX_ENABLED: 'true']
    }

    def 'should have default flake ref'() {
        given:
        def nix = new NixConfig([:], [:])
        expect:
        nix.flakeRef() == 'nixpkgs'
    }

    def 'should check flake ref'() {
        given:
        def nix = new NixConfig([flakeRef: 'github:NixOS/nixpkgs/nixos-24.05'], [:])
        expect:
        nix.flakeRef() == 'github:NixOS/nixpkgs/nixos-24.05'
    }

    def 'should check install options'() {
        given:
        def nix = new NixConfig([installOptions: '--impure'], [:])
        expect:
        nix.installOptions() == '--impure'
    }

    def 'should have default create timeout'() {
        given:
        def nix = new NixConfig([:], [:])
        expect:
        nix.createTimeout().toMillis() == 20 * 60 * 1000
    }
}
