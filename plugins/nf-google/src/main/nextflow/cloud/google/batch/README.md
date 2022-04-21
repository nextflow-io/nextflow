# Google Batch

### Get started 

0. Setup creds 

```
unset GOOGLE_APPLICATION_CREDENTIAL
unset GOOGLE_ACCESS_TOKEN
export SUPPRESS_GCLOUD_CREDS_WARNING=true
gcloud auth application-default login
```


1. Create this config file 

```
cat <<EOT > google.config
params.transcriptome = 'gs://rnaseq-nf/data/ggal/transcript.fa'
params.reads = 'gs://rnaseq-nf/data/ggal/gut_{1,2}.fq'
params.multiqc = 'gs://rnaseq-nf/multiqc'
params.outdir = 'gs://rnaseq-nf/results/rnaseq'

process.executor = 'google-batch'
process.container = 'quay.io/nextflow/rnaseq-nf:v1.1'
workDir = 'gs://rnaseq-nf/scratch'

google.project = 'my-nf-project-261815'
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

cat <<EOT > run.json
{
     "taskGroups": [
         {
             "taskSpec": {
                 "runnables": [
                     {
                         "container": {
                             "commands": [ "echo", "Hello" ],
                             "imageUri": "bash"
                         }
                     }
                 ],
                 "computeResource": {
                     "cpuMilli": 1000
                 }
             }
         }
     ],
     "logsPolicy": {
         "destination": "CLOUD_LOGGING"
     }
}
EOT

» gcurl --data @run.json https://${BatchAPI}/projects/${ProjectID}/locations/${Location}/jobs?job_id=c10


#### Plain tasks 

curl --header "Content-Type: application/json" \
    --header "Authorization: Bearer $GOOGLE_ACCESS_TOKEN" \
    https://batch.googleapis.com/v1alpha1/projects/${ProjectID}/locations/${Location}/jobs \
    | jq .jobs[].name -r


curl --header "Content-Type: application/json" \
    --header "Authorization: Bearer $GOOGLE_ACCESS_TOKEN" \
    https://batch.googleapis.com/v1alpha1/projects/${ProjectID}/locations/${Location}/jobs?job_id=SOME_TASK_ID \
    --data @run.json

  
#### Track log files

» CLOUDSDK_PYTHON_SITEPACKAGES=1 gcloud beta logging tail projects/my-nf-project-261815/logs/batch_task_logs
   

### Misc 

* Logging: https://github.com/googleapis/java-logging
