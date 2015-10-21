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

package nextflow.scm

/**
 * Local storage asset provider
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalRepositoryProvider extends RepositoryProvider {

    LocalRepositoryProvider( String project, ProviderConfig config) {
        assert project
        assert config
        assert config.path

        this.project = project
        this.config = config
        this.path = new File(config.path)
    }

    private File path

    @Override
    String getName() { 'local' }

    @Override
    String getEndpointUrl() {
        return "file:${new File(path, project)}"
    }

    @Override
    String getContentUrl(String path) {
        def base = new File(this.path, project)
        def result = new File(base, path)
        return "file:$result"
    }

    @Override
    String getCloneUrl() {
        return "file:${new File(path, project)}/.git"
    }

    @Override
    String getRepositoryUrl() {
        new File(path, project).toString()
    }

    @Override
    protected byte[] readBytes(String path) {
        return new File(new File(this.path, project), path).readBytes()
    }
}
