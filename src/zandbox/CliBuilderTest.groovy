
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

def args = ['-debug', 'xxx,yyy', '-t','-a','-l','*.groovy','--ciao','--bua','1']

def cli = new CliBuilder(usage:'ls')
cli.a('display all files')
cli.l('use a long listing format')
cli.t(longOpt:null, 'sort by modification time')
cli.debug('define the packages to debug')
def options = cli.parse(args)

assert options // would be null (false) on failure
assert options.a && options.l && options.t
assert options.debug
assert  options.arguments() == ['*.groovy','--ciao','--bua','1']