FROM openjdk:8-jre-alpine
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

# Add required dependencies
# - bash 
# - gnu coreutils 
# - curl 
RUN apk update && apk add bash && apk add coreutils && apk add curl && apk add graphviz

# see https://blogs.oracle.com/java-platform-group/java-se-support-for-docker-cpu-and-memory-limits
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow

# copy docker client
COPY dist/docker /usr/local/bin/docker
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
