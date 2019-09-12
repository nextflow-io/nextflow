build:  
	docker build -t nextflow/nextflow_dev -f Dockerfile ../

run:
	docker run --rm --name nextflow_dev -d -v $(PWD)/../:/nextflow_dev nextflow/nextflow_dev sleep infinity

exec:
	docker exec -it nextflow_dev bash

kill:
	docker kill nextflow_dev