FROM ubuntu

RUN apt-get update && apt-get -y install curl unzip && apt-get clean


RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && ./aws/install && rm -rf aws*

ADD fake_aws /fake_aws

ENV PATH=/fake_aws/bin/:$PATH
