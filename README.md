# ICGC-TCGA-PanCancer BWA Workflow

## Overview

This is the workflow for the TCGA/ICGC PanCancer project that aligns
whole genome sequences with BWA-Mem.

For more information about the project overall see the
[PanCancer wiki space](https://wiki.oicr.on.ca/display/PANCANCER/PANCANCER+Home).

More detailed documentation about the production use of this workflow can be
found in the [PanCancer-Info](https://github.com/ICGC-TCGA-PanCancer/pancancer-info)
project where we maintain our production documentation and SOPs.

## Building the Worklfow Docker Image

You can also build a Docker image that has the workflow ready to run in it.

    docker build -t pancancer/pcawg-bwa-workflow:2.6.7 .


## Running the Workflow with Cromwell

    cromwell run bwa-workflow.wdl bwa-workflow.json


## Sample Data

Some synthetic sample data.

* https://s3.amazonaws.com/oicr.workflow.bundles/released-bundles/synthetic_bam_for_GNOS_upload/hg19.chr22.5x.normal2.bam
* https://s3.amazonaws.com/oicr.workflow.bundles/released-bundles/synthetic_bam_for_GNOS_upload/hg19.chr22.5x.normal.bam

## Reference Data

We use a specific reference based on GRCh37.

* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.fai
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.64.amb
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.64.ann
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.64.bwt
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.64.pac
* http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.64.sa

## Authors

* Brian O'Connor <boconnor@oicr.on.ca>
* Junjun Zhang <Junjun.Zhang@oicr.on.ca>
* Adam Wright <adam.wright@oicr.on.ca>

## Contributors

* Keiran Raine: PCAP-Core and BWA-Mem workflow design
* Roshaan Tahir: Original BWA-Align workflow design
* Adam Struck: WDL implementation
