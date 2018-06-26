#!/bin/bash
# TODO: This is a temporary workaround to deploy local nextflow binary go GCE bucket

gsutil cp build/releases/nextflow-*-all gs://nextflow-nextcode-dev/
