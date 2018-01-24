FROM openjdk:8-jre-alpine
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

# Add required dependencies
# - bash 
# - gnu coreutils 
# - curl 
RUN apk update && apk add bash && apk add coreutils && apk add curl

# download runtime
RUN curl https://www.nextflow.io/releases/v0.27.1/nextflow > /usr/local/bin/nextflow \
  && chmod 755 /usr/local/bin/nextflow \
  && sync \
  && nextflow -download

# see https://blogs.oracle.com/java-platform-group/java-se-support-for-docker-cpu-and-memory-limits
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap'

# define the entry point
ENTRYPOINT ["nextflow"] 
