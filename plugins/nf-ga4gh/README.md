# GA4GH plugin for Nextflow

This plugin implements the support for GA4GH APIs for Nextflow. Currently only supports the [Task Execution Service (TES) API](https://github.com/ga4gh/task-execution-schemas).

## SDK Generation

[Swagger Codegen](https://github.com/swagger-api/swagger-codegen) was used to generate the Java SDK for TES based on the [TES OpenAPI specification](https://github.com/ga4gh/task-execution-schemas/blob/develop/openapi/task_execution_service.openapi.yaml).

The easiest way to generate the Java SDK is with the [Docker image](https://github.com/swagger-api/swagger-codegen#swagger-codegen-cli-docker-image):

```bash
# download the TES OpenAPI spec
wget https://github.com/ga4gh/task-execution-schemas/raw/v1.1/openapi/task_execution_service.openapi.yaml

# convert the spec to JSON
python3 -c 'import sys, yaml, json; y=yaml.safe_load(sys.stdin.read()); print(json.dumps(y, indent=2, default=str))' \
  < task_execution_service.openapi.yaml \
  > task_execution_service.openapi.json

# generate the Java SDK
docker run -v ${PWD}:/local swaggerapi/swagger-codegen-cli-v3 generate -i /local/task_execution_service.openapi.json -l java -o /local
```
