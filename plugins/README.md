# Nextflow plugins

This directory should contain plugin subprojects for Nextflow. 

The `build.gradle` defines the main actions for each plugin.

## Plugin subproject structure 
 
### Plugin structure 

The plugin subproject defines its own `build.gradle` and setup the required dependencies. 

Minimal dependencies shown below: 

``` 
dependencies {
    compileOnly project(':nextflow')
    compileOnly 'org.slf4j:slf4j-api:1.7.10'
    compileOnly 'org.pf4j:pf4j:3.4.1'

    testImplementation project(':nextflow')
    testImplementation "org.codehaus.groovy:groovy:3.0.5"
    testImplementation "org.codehaus.groovy:groovy-nio:3.0.5"
}
``` 

Each plugin subproject directory name has to begin with the prefix `nf-` and must include 
a file named `src/resources/META-INF/MANIFEST.MF` which contains the plugin metadata. 
The manifest content looks like the following:

```
Manifest-Version: 1.0
Plugin-Class: the.plugin.ClassName
Plugin-Id: the-plugin-id
Plugin-Provider: Some Provider Name
Plugin-Version: 0.0.0
```   
  
## Environment variables 

* `NXF_PLUGINS_MODE`: Define the plugin system execution mode, either *prod* for production or *dev* for development
  (see below for details).  
* `NXF_PLUGINS_DIR`: the path where the plugins archives are stored/loaded. Default is `$NXF_HOME/plugins` for 
  production mode and `$PWD/plugins` for dev mode. 
* `NXF_PLUGINS_DEFAULT`: Whenever use the default plugins when no plugin is specified in the config file.   
* `NXF_PLUGINS_DEV: Comma separate separated list of development plugins root directories

## Development environment

When running in development the plugin system uses the `DevPluginClasspath` to load plugins classes 
from each plugin project build path e.g. `$PWD/plugins/nf-amazon/build/classes` and 
`$PWD/plugins/nf-amazon/build/target/libs` (for deps libraries).    


## The plugins repository 

The plugins meta-info are published via a GitHub repository at [https://github.com/nextflow-io/plugins]()
and accessible through the URL [https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json]().

The repository index has the following structure: 

```
[
  {
    "id": "nf-amazon",
    "releases": [
      {
        "version": "0.2.0",
        "url": "https://github.com/nextflow-io/nf-amazon/releases/download/0.2.0/nf-amazon-0.2.0.zip",
        "date": "2020-10-12T10:05:44.28+02:00",
        "sha512sum": "9e9e33695c1a7c051271..."
      }
    ]
  },
  :
]
```     


## Plugins 

A plugin is a ZIP file holding either the plugin classes and the required dependencies JAR file.  

Nextflow core plugins are stored in the corresponding GitHub project release page. However, this is not a
strict requirement, it has been chosen to simplify the build deployment process and provide a more consistent 
download experience keeping all of them with the GitHub [nextflow-io](https://github.com/nextflow-io) organization.   

### The installation process 

Plugins need to be declared in the `nextflow.config` file using the plugins scope, eg. 

```
plugins {
    id 'nf-amazon@0.2.0'
}
```     

If the plugins is not locally available Nextflow check in the repository index for the download URL, 
download the ZIP in a temporary file, unzip and store the final plugin in the directory specified by the 
variable `NXF_PLUGINS_DIR` (default: `$NXF_HOME/plugins`). 

Finally, since each Nextflow run can have a different set of plugins (and version) requirement, each Nextflow 
instance keep local plugins directory root in directory `$PWD/.nextflow/plr/<unique id>` symlinking the exact list 
of plugins directory required for the current Nextflow instance.

If no plugins are specified in the nextflow.config file, Nextflow default plugins are automatically added. 
The default plugins list is defined in the Nextflow resources file included in the distribution runtime 
`./modules/nextflow/src/main/resources/META-INF/plugins-info.txt`. 

To disable the use of defualt plugins set the following variable `NXF_PLUGINS_DEFAULT=false`.

## Gradle Tasks 

### makeZip
    
Creates the plugin the zip file and the json meta file in the
subproject `build/libs` directory.

```
» ls -l1 $PWD/plugins/nf-tower/build/libs/
nf-tower-0.1.0.jar
nf-tower-0.1.0.json
nf-tower-0.1.0.zip
```               

### copyPluginLibs

Copies plugin dependencies jar files in the plugin build directory ie. `$PWD/plugins/nf-amazon/build/target/libs`. 
This is only needed when launching the plugin in *development* mode. 

### copyPluginZip

Copies the plugin zip file to the root project build dir ie. `$PWD/build/plugins/`.

### uploadPlugin

Uploads the plugin ZIP and meta (JSON) files to the corresponding GitHub repository. Options: 

* `release`: the plugin version e.g. `1.0.1`
* `repo`: the GitHub repository name e.g. `nf-amazon`
* `owner`: the GitHub owning organization e.g. `nextflow-io`
* `skipExisting`: do not upload a file already existing (only if the checksum is the same, default: `true`). 
* `dryRun`: execute the tasks without uploading file (default: `false`).
* `overwrite`: prevent to overwrite a remote file already existing (default: `false`).
* `userName`: the user name for authenticate GitHub API requests
* `authToken`: the personal token to authenticate GitHub API requests  

### upload

Uploads the plugin both zip and meat files. 

### publishIndex

Upload the plugins index to the repository hosted at [https://github.com/nextflow-io/plugins](https://github.com/nextflow-io/plugins), which makes 
them accessible through the URL [https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json](https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json). 


## Links 

* https://pf4j.org/
* https://proandroiddev.com/understanding-gradle-the-build-lifecycle-5118c1da613f
