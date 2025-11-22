FROM amazoncorretto:21-al2023
RUN yum install -y procps-ng shadow-utils which

ENV NXF_HOME=/.nextflow
ARG TARGETPLATFORM=linux/amd64

# copy docker client
COPY dist/${TARGETPLATFORM}/docker /usr/local/bin/docker
COPY entry.sh /usr/local/bin/entry.sh
COPY nextflow /usr/local/bin/nextflow

# download runtime
RUN mkdir /.nextflow \
 && chmod 755 /usr/local/bin/nextflow \
 && chmod 755 /usr/local/bin/entry.sh \
 && nextflow info

# define the entry point
ENTRYPOINT ["/usr/local/bin/entry.sh"]
