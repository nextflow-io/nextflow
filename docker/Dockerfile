FROM groovy:jdk8 as build

USER root
RUN apt-get update \
    && apt-get install -yq make
#USER groovy

WORKDIR /src
COPY . /src

RUN make deps
RUN make compile

#RUN make test
#RUN make smoke

RUN make pack

#RUN make install

RUN cd docker && make dist/docker

FROM openjdk:8-jre-alpine
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

# Add required dependencies
# - bash 
# - gnu coreutils 
# - curl 
RUN apk --no-cache add bash coreutils curl

# see https://blogs.oracle.com/java-platform-group/java-se-support-for-docker-cpu-and-memory-limits
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow

# copy docker client
COPY --from=build /src/docker/dist/docker /usr/local/bin/docker
COPY docker/entry.sh /usr/local/bin/entry.sh
COPY --from=build /src/nextflow /usr/local/bin/nextflow

# download runtime
RUN mkdir /.nextflow \
 && touch /.nextflow/dockerized \
 && chmod 755 /usr/local/bin/nextflow \
 && chmod 755 /usr/local/bin/entry.sh \
 && nextflow info

# define the entry point
ENTRYPOINT ["/usr/local/bin/entry.sh"]
