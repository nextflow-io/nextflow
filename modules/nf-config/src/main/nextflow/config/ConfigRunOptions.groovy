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
package nextflow.config

/**
 * Interface to get `run` command options used for configuration purposes.
 *
 * @author Jorge Ejarque (jorge.ejarque@seqera.io)
 */
interface ConfigRunOptions {
    String getBucketDir()
    String getCloudCachePath()
    String getLibPath()
    String getOutputDir()
    String getWorkDir()
    Boolean getCacheable()
    boolean getPreview()
    String getResume()
    String getRunName()
    String getTest()
    boolean getStubRun()
    Integer getPoolSize()
    long getPollInterval()
    Integer getQueueSize()
    String getPlugins()
    String getProfile()
    List<String> getRunConfig()
    Map<String,String> getExecutorOptions()
    Map<String,String> getProcess()
    Map<String,String> getClusterOptions()
    Map<String,String> getEnv()
    boolean getExportSysEnv()
    ConfigCliOptions getConfigCliOptions()
    Boolean getWithoutConda()
    String getWithConda()
    Boolean getWithoutSpack()
    String getWithSpack()
    String getWithTrace()
    String getWithReport()
    String getWithTimeline()
    String getWithDag()
    String getWithNotification()
    String getWithWebLog()
    String getWithTower()
    String getWithWave()
    String getWithFusion()
    boolean getWithoutDocker()
    def getWithDocker()
    def getWithPodman()
    def getWithSingularity()
    def getWithApptainer()
    def getWithCharliecloud()
    String getDumpHashes()
    String getDumpChannels()
    String getNormalizeResumeId(String id)

}