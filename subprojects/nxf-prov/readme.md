# Nextflow implements Research Object specification
## Google Summer of Code 2018 Project - OBF

With this new subproject, we want to allow Nextflow users to generate a zip file of their pipeline following the specifications of a Research Object.

# What is a Research Object?
Research Object is a method for the identification, aggregation, and exchange of information. The primary goal is to provide a way to associate together related resources from the same project, like the pipeline, auxiliar scripts, data, slides or the final article.
The Research Object idea is motivated by a desire to improve reproducibility. The main three principles of it are: 
*Identity*, providing a unique identifier to the project, as the DOI for the publication or the ORCID for the scientist. *Aggregation*, allowing the author to wrap all the elements used for the project, slides, the article, scripts … With the Research Object, we can share all the elements of a project together with the same ID. And finally, the last main principle of Research Object is the *Annotation*, it means to provide an extra metadata to know the relation between elements, when and how they were produced.
# How to use it
First of all, you will need to follow a couple of steps, to allow Nextflow to generate the RO with the information it needs.
You will need to include the directory of the output files. Inside the ```nextflow.config``` file, you will need to add/include the author's name (optional), the author's ORCID (optional)and the output directory (mandatory).  
```sh
manifest {
  description = '...'
  author = 'autorFOO'
  ORCID = 'orcidFOO'
  outdir = './results'
}
```
This procedure is under development to improve it and get the information in an automatic way.

Once you have completed the information needed by Nextfow, you need to have your pipeline executed. Usually, it will be something like this: ```nextflow run <your main file> --parameters <your pipeline parameters>```
Once the pipeline is done. You will be able to generate the RO zip file with the command line ```./prov.sh```
Then, it will appear on the console some information about the procedure, and it will show some *WARN* if we missed some information to add or there is any unexpected problem during the execution.
 ## How it works
 The way it works it's easy. During the execution of the pipeline, Nextflow is capturing all the metadata and information on the fly, and it will store it on the cache and on a provenance file.
 Once the user runs the ```./prov.sh```, Nextflow starts converting all the information and metadata into the RO object.
 First step, is to generate the file structure in the zip file. This file will contain the following folders:
 - **Data** 
 - **Workflow** 
 - **Snapshot** 
 - **Output** 
 - **Metadata** 
 - **.ro** 

![Structure v2.1](https://files.gitter.im/privateEdgano/Lobby/NnQB/image.png)

Some of the dummy examples are displaied bellow:

Snapshot-> **commandLine.txt** :
```
./launch.sh run nextflow-io/rnatoy -with-docker
```
**mimetype**:
```
application/vnd.wf4ever.robundle+zip
```
Metadata->**provenance.json**:
```
"wasGeneratedBy": {
    "PROV:generatedBy_2_/Users/edgargarriga/CBCRG/nextflow/work/e1/647ba7f1b1015ca8356c35cfe256bf/accepted_hits.bam": {
      "prov:entity": "PROV:/Users/edgargarriga/CBCRG/nextflow/work/e1/647ba7f1b1015ca8356c35cfe256bf/accepted_hits.bam",
      "prov:role": {
        "$": "",
        "type": "xsd:string"
      },
      "prov:activity": "PROV:activity_2"
    },
    ...
```
```
"activity": {
    "PROV:activity_3": {
      "prov:startTime": "2018-08-10T17:48:37.253+02:00",
      "prov:type": {
        "$": "Process",
        "type": "PROV:activityType"
      },
      "prov:label": "mapping (ggal_liver)",
      "prov:endTime": "2018-08-10T17:49:00.684+02:00"
    },
    ...
```
```
"wasAssociatedWith": {
    "PROV:associatedWith_1": {
      "prov:agent": "PROV:agent_activity_1",
      "prov:activity": "PROV:activity_1"
    },
    ...
```
```
"used": {
    "PROV:used_3_/Users/edgargarriga/.nextflow/assets/nextflow-io/rnatoy/data/ggal/ggal_liver_2.fq": {
      "prov:entity": "PROV:/Users/edgargarriga/.nextflow/assets/nextflow-io/rnatoy/data/ggal/ggal_liver_2.fq",
      "prov:activity": "PROV:activity_3"
    },
    ...
```
```
"entity": {
    "PROV:/nextflow/work/bb/2b5f718ecb9176cee3cf22f4dbdfcb/genome.index.4.bt2": [
      {
        "prov:value": {
          "$": "genome.index.4.bt2",
          "type": "PROV:fileName"
        },
        "prov:type": [
          {
            "$": "7667affc50c9a2f1f6cdbe4523b89a007f79388d8636a12058e5bbccade0a1aa",
            "type": "PROV:SHA256"
          },
          {
            "$": "42751",
            "type": "PROV:fileSize"
          }
        ]
      },
      ...
```
Metadata->**metadata.xml**:
```
Command Line: ./launch.sh run nextflow-io/rnatoy -with-docker
UUID: cd30a6ba-873b-485f-9773-ca1060d46241
Nextflow version: 0.31.0.4891
```
Metadata->**log.txt.**:
```
08/01/2018 12:03:58 -- autorFOO -- 399b7fc9-7a2f-4754-ba33-27450fc43588 -- ./launch.sh run nextflow-io/rnatoy -with-docker -resume
08/01/2018 12:04:14 -- autorFOO -- 399b7fc9-7a2f-4754-ba33-27450fc43588 -- ./launch.sh run nextflow-io/rnatoy -with-docker -resume
08/01/2018 12:04:26 -- autorFOO -- 399b7fc9-7a2f-4754-ba33-27450fc43588 -- ./launch.sh run nextflow-io/rnatoy -with-docker -resume
...
```
.ro->**manifest.json**:
```
{
  "@context" : [ "https://w3id.org/bundle/context" ],
  "id" : "/",
  "manifest" : [ "manifest.json" ],
  "createdOn" : "2018-08-10T15:50:53.769Z",
  "createdBy" : {
    "orcid" : "**ORCID_not_provided**",
    "name" : "autorFOO"
  },
  "aggregates" : [ {
    "uri" : "/outputs/",
    "createdOn" : "2018-08-10T15:50:53.768Z",
    "bundledAs" : {
      "uri" : "urn:uuid:51c1bf0e-980e-49e2-bb23-eee34536cbeb",
      "folder" : "/"
    }
  },
  ...
```
This will be the structure inside the zip file. Then each folder will contain a series of files with a determined information and a specific rol on the RO.
Inside the Data folder will be all the inputs files used by the pipeline. Then, inside the Workflow folder we are able to see all the files of the nextfow pipeline directory, with the NF script, the config file, the auxiliar scripts, etc...
On the snapshot we will find a sh file the the command line used to ran the pipeline. It will allow us to run the same NF version with the exact same parameters. On the Output folder we will find all the output files/results of our pipeline. Inside the Metadata folder, we can find some interesting files like the provenance.json. This files is capturing all the relations between channels and process. Being able to track wich process generated which file, when it was done and how long it took. We are able to see wich command line was used for each file and wich container technology (docker/singularity).
ON the metadata.xml file, we are able to see the exact command line used to run the pipeline, the UUID of the run (generated by nextflow) and finally we are able to know the version of nextflow used to run all the process. On the log.txt file we can see when and who generate older versions of the same RO zip file. 
Finally, inside the .ro folder, we can find the manifest.json file with all the information and metainformation of the zip file. When it was generated, which files are included, the uuid of the Research Object zip file.  