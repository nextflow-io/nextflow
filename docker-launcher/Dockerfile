FROM docker.io/pditommaso/nf-lancher:j17-base
COPY .nextflow /.nextflow
COPY nextflow /usr/bin/nextflow
ENV NXF_HOME=/.nextflow
RUN nextflow info
RUN NXF_PLUGINS_DEFAULT=false nextflow plugin install nf-tower,nf-wave,nf-cloudcache,nf-azure,nf-google,nf-amazon,xpack-amzn,xpack-google,nf-cloudcache
RUN rm /.nextflow/plugins/nf-azure-1.3.3-patch2/lib/gson-2.8.0.jar


