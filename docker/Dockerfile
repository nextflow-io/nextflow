FROM amazoncorretto:17.0.3
RUN yum install -y procps-ng

FROM amazoncorretto:17.0.3
COPY --from=0 /bin/ps /bin/ps
COPY --from=0 /lib64/libprocps.so.4 /lib64/libprocps.so.4
COPY --from=0 /usr/lib64/libsystemd.so.0 /usr/lib64/libsystemd.so.0
ENV NXF_HOME=/.nextflow
ARG TARGETPLATFORM

# copy docker client
COPY dist/${TARGETPLATFORM}/docker /usr/local/bin/docker
COPY entry.sh /usr/local/bin/entry.sh
COPY nextflow /usr/local/bin/nextflow

# download runtime
RUN mkdir /.nextflow \
 && touch /.nextflow/dockerized \
 && chmod 755 /usr/local/bin/nextflow \
 && chmod 755 /usr/local/bin/entry.sh \
 && nextflow info

# define the entry point
ENTRYPOINT ["/usr/local/bin/entry.sh"]
