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

package io.seqera.tower.plugin.dataset

import io.seqera.tower.plugin.exception.ForbiddenException
import io.seqera.tower.plugin.exception.NotFoundException
import io.seqera.tower.plugin.exception.UnauthorizedException

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException
import java.time.OffsetDateTime

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.model.OrgAndWorkspaceDto
import io.seqera.tower.plugin.TowerClient
import nextflow.exception.AbortOperationException

/**
 * Typed client for Seqera Platform dataset API endpoints.
 * Delegates HTTP execution to {@link TowerClient#sendApiRequest}, inheriting its
 * authentication token management and retry policy without exposing the underlying
 * HTTP client.
 *
 * @author Seqera Labs
 */
@Slf4j
@CompileStatic
class SeqeraDatasetClient {

    private final TowerClient towerClient

    SeqeraDatasetClient(TowerClient towerClient) {
        this.towerClient = towerClient
    }

    private String getEndpoint() { towerClient.endpoint }

    /**
     * @return current user info (id, userName, etc.) from GET /user-info
     */
    Long getUserId() {
        try {
             final info = towerClient.getUserInfo()
             if( info?.id == null )
                 throw new AbortOperationException("Unable to retrieve user ID from Seqera Platform — check your access token")
             return info.id as long
        }catch( UnauthorizedException e ){
            throw new AbortOperationException(e.getMessage())
        }catch( ForbiddenException e){
            throw new AccessDeniedException("${endpoint}/user-info", null, e.message)
        }catch(NotFoundException e){
            throw new NoSuchFileException("${endpoint}/user-info")
        }
    }

    /**
     * @return all orgs and workspaces accessible to the given user from GET /user/{userId}/workspaces
     */
    List<OrgAndWorkspaceDto> listUserWorkspacesAndOrgs(long userId) {
        try {
            final list = towerClient.listUserWorkspacesAndOrgs(userId as String)
            return list.collect { m -> mapOrgAndWorkspace(m) }
        } catch( UnauthorizedException e ){
            throw new AbortOperationException(e.getMessage())
        } catch( ForbiddenException e){
            throw new AccessDeniedException("${endpoint}/user/$userId/workspaces", null, e.message)
        } catch(NotFoundException e){
            throw new NoSuchFileException("${endpoint}/user/$userId/workspaces")
        }
    }



    /**
     * @return all datasets in the given workspace from GET /datasets?workspaceId={workspaceId}
     */
    List<DatasetDto> listDatasets(long workspaceId) {
        final url = "${endpoint}/datasets?workspaceId=${workspaceId}"
        log.debug "SeqeraDatasetClient GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final list = json.datasets as List<Map>
        return list ? list.collect { m -> mapDataset(m) } : Collections.<DatasetDto>emptyList()
    }

    /**
     * Create a new dataset in the given workspace via POST /datasets?workspaceId={workspaceId}.
     * @return the created dataset DTO
     */
    DatasetDto createDataset(long workspaceId, String name) {
        final url = "${endpoint}/datasets?workspaceId=${workspaceId}"
        log.debug "SeqeraDatasetClient POST $url name=$name"
        final resp = towerClient.sendApiRequest(url, [name: name], 'POST')
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        return mapDataset(json.dataset as Map)
    }

    /**
     * @return all versions for the given dataset from GET /datasets/{datasetId}/versions
     */
    List<DatasetVersionDto> listVersions(String datasetId, long workspaceId) {
        final url = "${endpoint}/datasets/${datasetId}/versions?workspaceId=${workspaceId}"
        log.debug "SeqeraDatasetClient GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final list = json.versions as List<Map>
        return list ? list.collect { m -> mapVersion(m) } : Collections.<DatasetVersionDto>emptyList()
    }

    /**
     * Download a dataset version as an InputStream.
     * GET /datasets/{datasetId}/v/{version}/n/{fileName}
     * The fileName must exactly match DatasetVersionDto.fileName from upload time.
     */
    InputStream downloadDataset(String datasetId, String version, String fileName, long workspaceId) {
        final encodedName = new URI(null, null, fileName, null).rawPath
        final url = "${endpoint}/datasets/${datasetId}/v/${version}/n/${encodedName}?workspaceId=$workspaceId"
        log.debug "SeqeraDatasetClient GET $url (streaming)"
        try {
            return towerClient.sendStreamingRequest(url)
        }
        catch (UnauthorizedException e) {
            throw new AbortOperationException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        }
        catch (ForbiddenException e) {
            throw new AccessDeniedException(url, null, "Forbidden — check workspace permissions")
        }
        catch (NotFoundException e) {
            throw new NoSuchFileException(url)
        }
    }

    // ---- private helpers ----

    private static void checkFsResponse(TowerClient.Response resp, String url) {
        if (!resp.error) return
        final code = resp.code
        if (code == 401)
            throw new AbortOperationException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        if (code == 403)
            throw new AccessDeniedException(url, null, "Forbidden — check workspace permissions")
        if (code == 404)
            throw new NoSuchFileException(url)
        throw new IOException("Seqera API error: HTTP ${code} for ${url}")
    }

    private static OrgAndWorkspaceDto mapOrgAndWorkspace(Map m) {
        final dto = new OrgAndWorkspaceDto()
        dto.orgId = (m.orgId as Long) ?: 0L
        dto.orgName = m.orgName as String
        dto.workspaceId = (m.workspaceId as Long) ?: 0L
        dto.workspaceName = m.workspaceName as String
        dto.workspaceFullName = m.workspaceFullName as String
        return dto
    }

    private static DatasetDto mapDataset(Map m) {
        final dto = new DatasetDto()
        dto.id = m.id as String
        dto.name = m.name as String
        dto.description = m.description as String
        dto.version = (m.version as Long) ?: 0L
        dto.mediaType = m.mediaType as String
        dto.workspaceId = (m.workspaceId as Long) ?: 0L
        dto.dateCreated = m.dateCreated ? OffsetDateTime.parse(m.dateCreated as String) : null
        dto.lastUpdated = m.lastUpdated ? OffsetDateTime.parse(m.lastUpdated as String) : null
        return dto
    }

    private static DatasetVersionDto mapVersion(Map m) {
        final dto = new DatasetVersionDto()
        dto.datasetId = m.datasetId as String
        dto.version = (m.version as Long) ?: 0L
        dto.fileName = m.fileName as String
        dto.mediaType = m.mediaType as String
        dto.hasHeader = (m.hasHeader as Boolean) ?: false
        dto.dateCreated = m.dateCreated ? OffsetDateTime.parse(m.dateCreated as String) : null
        dto.disabled = (m.disabled as Boolean) ?: false
        return dto
    }
}
