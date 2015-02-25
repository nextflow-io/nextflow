
/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

def conf1 = '''
    env {
        alpha = 'env1'
        beta = 'env1'
        HOME = "$HOME:xxx"
    }

    '''

def conf2 = '''
    task {
        processor = 'local'
        echo = true
        shell = 'bash'
        validExitCodes = 0
        maxForks = 1
        maxDuration = '10 m'
        maxMemory = '1G'
        maxInactive = '1 m'
        errorStrategy = 'xx'
    }

    env {
        beta = 'env2'
        delta = 'env2'
        HOME = "$HOME:yyy"
    }

    x = "$HOME"

    '''

def builder = new StringBuilder("env {\n")
System.getenv().sort().each { name, value -> builder << "${name.replace('.','_')}='$value'\n" }
builder << "}"

def slurp0 = new ConfigSlurper()
slurp0.setBinding(System.getenv())

def slurp1 = new ConfigSlurper()
slurp1.setBinding(System.getenv())

def slurp2 = new ConfigSlurper()
slurp2.setBinding(System.getenv())


def config0 = slurp0.parse(builder.toString())
def config1 = slurp1.parse(conf1)
def config2 = slurp2.parse(conf2)

config0.merge(config1)
config0.merge(config2)


assert config0.env.alpha == 'env1'
assert config0.env.beta == 'env2'
assert config0.env.delta == 'env2'
assert config0.env.HOME == System.getenv()['HOME'] + ':yyy'
assert config0.x == System.getenv()['HOME']


