FROM openjdk:8-jre-alpine
MAINTAINER Paolo Di Tommmaso <paolo.ditommaso@gmail.com> 

# Add requiders dependencies 
# - bash 
# - gnu coreutils 
# - curl 
RUN apk update && apk add bash && apk add coreutils && apk add curl

# pull deps 
RUN curl https://www.nextflow.io/releases/v0.25.2/nextflow > /usr/local/bin/nextflow \
  && chmod 755 /usr/local/bin/nextflow \
  && sync \
  && nextflow -download

# define the entry point 
ENTRYPOINT nextflow 
