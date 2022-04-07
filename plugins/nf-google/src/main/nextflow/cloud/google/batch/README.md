# Google Batch


* Logging: https://github.com/googleapis/java-logging


### Usage 


export ProjectID=my-nf-project-261815
export BatchAPI=batch.googleapis.com/v1alpha1
export Location=us-central1
alias gcurl='curl --header "Content-Type: application/json" --header "Authorization: Bearer $(gcloud auth print-access-token)"'


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

» gcurl https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs/container-425
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
