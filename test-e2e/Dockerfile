FROM amazoncorretto:21-al2023
COPY .nextflow /.nextflow
COPY nextflow /usr/local/bin/nextflow
ENV NXF_HOME=/.nextflow
RUN nextflow info
RUN NXF_PLUGINS_DEFAULT=false nextflow plugin install xpack-amzn,xpack-google

FROM scratch
COPY --from=0 /.nextflow /.nextflow
COPY --from=0 /usr/local/bin/nextflow /usr/local/bin/nextflow
