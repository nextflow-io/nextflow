# Google Batch

### Get started 

0. Setup access token 

```
export GOOGLE_ACCESS_TOKEN=$(gcloud auth print-access-token)
```

1. Create this config file 

```
cat <<EOT > google.config
params.transcriptome = 'gs://rnaseq-nf/data/ggal/transcript.fa'
params.reads = 'gs://rnaseq-nf/data/ggal/gut_{1,2}.fq'
params.multiqc = 'gs://rnaseq-nf/multiqc'
process.executor = 'google-batch'
workDir = 'gs://rnaseq-nf/scratch'
google.region  = 'us-central1'
EOT
```
    
2. Compile 

```
make compile
```

3. Run it 

```
./launch.sh run pditommaso/hello -c google.config 
```

### Config 

export ProjectID=my-nf-project-261815
export BatchAPI=batch.googleapis.com/v1alpha1
export Location=us-central1
export GOOGLE_ACCESS_TOKEN=$(gcloud auth print-access-token)
alias gcurl='curl --header "Content-Type: application/json" --header "Authorization: Bearer $GOOGLE_ACCESS_TOKEN"'

#### List jobs 

» gcurl https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs | jq .jobs[].name -r
  
#### Delete a job 

» gcurl -X DELETE https://${BatchAPI}/$JOB_NAME

#### Delete all jobs

for x in $(gcurl  https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs | jq .jobs[].name -r); do 
  echo "Delete job $x"
  gcurl -X DELETE https://${BatchAPI}/$x
  sleep 0.5
done

#### Submit a job 

» gcurl --data @run.json https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs?job_id=c10
{
  "name": "projects/my-nf-project-261815/locations/us-central1/jobs/container06",
  "uid": "j-23a3f19a-d06a-4206-9c38-62e0a237da7b",
  "taskGroups": [
    {
      "name": "projects/770252724889/locations/us-central1/jobs/container06/taskGroups/group0",
      "taskSpec": {
        "computeResource": {
          "cpuMilli": "1000",
          "memoryMib": "2048"
        },
        "volumes": [
          {
            "gcs": {
              "remotePath": "nf_bucket"
            },
            "mountPath": "/mnt/share"
          }
        ],
        "runnables": [
          {
            "container": {
              "imageUri": "busybox",
              "commands": [
                "-c",
                "echo \"task ${BATCH_TASK_INDEX}\" \u003e\u003e /mnt/share/batch-log${BATCH_TASK_INDEX}.txt"
              ],
              "entrypoint": "/bin/sh"
            }
          }
        ]
      },
      "taskCount": "1",
      "parallelism": "1",
      "taskCountPerNode": "1"
    }
  ],
  "status": {
    "state": "QUEUED"
  },
  "createTime": "2021-11-02T22:53:56.081702894Z",
  "updateTime": "2021-11-02T22:53:56.081702894Z",
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}

» gcurl --data @run.json https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs?job_id=c10
{
  "name": "projects/my-nf-project-261815/locations/us-central1/jobs/container06",
  "uid": "j-23a3f19a-d06a-4206-9c38-62e0a237da7b",
  "taskGroups": [
    {
      "name": "projects/770252724889/locations/us-central1/jobs/container06/taskGroups/group0",
      "taskSpec": {
        "computeResource": {
          "cpuMilli": "1000",
          "memoryMib": "2048"
        },
        "volumes": [
          {
            "gcs": {
              "remotePath": "nf_bucket"
            },
            "mountPath": "/mnt/share"
          }
        ],
        "runnables": [
          {
            "container": {
              "imageUri": "busybox",
              "commands": [
                "-c",
                "echo \"task ${BATCH_TASK_INDEX}\" \u003e\u003e /mnt/share/batch-log${BATCH_TASK_INDEX}.txt"
              ],
              "entrypoint": "/bin/sh"
            }
          }
        ]
      },
      "taskCount": "1",
      "parallelism": "1",
      "taskCountPerNode": "1"
    }
  ],
  "status": {
    "state": "FAILED",
    "taskGroups": {
      "group0": {
        "counts": {
          "FAILED": "1"
        }
      }
    }
  },
  "createTime": "2021-11-02T22:53:56.081702894Z",
  "updateTime": "2021-11-02T22:53:56.132247867Z",
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}


» CLOUDSDK_PYTHON_SITEPACKAGES=1 gcloud beta logging tail projects/my-nf-project-261815/logs/batch_task_logs
   

### Misc 

* Logging: https://github.com/googleapis/java-logging
