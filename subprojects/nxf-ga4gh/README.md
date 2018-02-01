# GA4GH API 

http://ga4gh.org/

## Task Execution Schema (TES)

### Schema definition 

https://github.com/ga4gh/task-execution-schemas/


## Workflow Execution Schema (WES)

### Schema definition 

https://github.com/ga4gh/workflow-execution-schemas


### Reference implementation 

https://github.com/common-workflow-language/workflow-service/tree/ga4gh-wes

### Funnel backend

 - Tested with release [0.2.0](https://github.com/ohsu-comp-bio/funnel/releases/tag/0.2.0)

 Funnel configuration file:

 ```
 Storage:
   Local:
     AllowedDirs:
       - /Users
 ```

 Funnel launch command line:

 ```
 ./funnel-darwin-amd64 server -c funnel.yml --http-port 8000
 ```

### Other 

* [Swagger editor](http://editor.swagger.io) 