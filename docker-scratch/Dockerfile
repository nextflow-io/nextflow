FROM amazoncorretto:17-al2023
COPY .nextflow /.nextflow
COPY nextflow /usr/bin/nextflow
ENV NXF_HOME=/.nextflow
RUN nextflow info
RUN NXF_PLUGINS_DEFAULT=false nextflow plugin install nf-tower,nf-wave,nf-cloudcache,nf-azure,nf-google,nf-amazon,xpack-amzn,xpack-google,nf-cloudcache

#FROM scratch
#COPY --from=0 /.nextflow /.nextflow
#COPY --from=0 /usr/bin/nextflow /usr/bin/nextflow

