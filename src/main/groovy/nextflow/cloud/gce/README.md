Google Compute Engine Cloud Driver
==================================

Usage example
-------------

Download google cloud CLI and initialize using
`gcloud init`

Create and download credentials for google cloud in the console
Go to APIs & Services
Credentials -> Create Credentials
Select Service account key.
Download json file and save as creds.json
Export as env variable:
```bash
export GOOGLE_APPLICATION_CREDENTIALS=$PWD/creds.json
```

Create nextflow.config file in project root directory. Example:
```
cloud {
  imageId = 'centos-cloud/global/images/centos-7-v20180815'
  instanceType = 'n1-standard-1'
}

gce {
  project = 'your-project-id'
  zone = 'us-central1-f'
}

```

Build project:
```bash
make compile pack install
```

Launch cluster:
```bash
./nextflow cloud create my-cluster -c 2 -driver gce 
```

This should launch master node and one slave

Google Pipelines API
====================

Usage Example
-------------
Download google cloud CLI and initialize using
`gcloud init`

Create and download credentials for google cloud in the console
Go to APIs & Services
Credentials -> Create Credentials
Select Service account key.
Download json file and save as creds.json
Export as env variable:
```bash
export GOOGLE_APPLICATION_CREDENTIALS=$PWD/creds.json
```

Create nextflow.config file in project root directory. The config must specify the following parameters:

* workDir - *Must be a GS bucket*
* process.executor - **googlepipelines**
* cloud.instanceType - *valid GCE instance type*
* gce.project - *gce project to run in*
* gce.zone - *gce zone to run in*

Example:
```
workDir = 'gs://<your bucket>/<directory>/'

process {
  executor = 'googlepipelines'
}

cloud {
  instanceType = 'n1-standard-1'
}

gce {
  project = 'your-project-id'
  zone = 'us-central1-f'
}

```

Note that all tasks defined in your nextflow scripts must define a container to run in, or you can configure it
globally with the *process.container* option in your nextflow.config

```
process.container = 'ubuntu:latest'
```


Build project:
```bash
make compile pack install
```

Run your nextflow script.


