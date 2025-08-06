#!/bin/bash

cap md5 --slurm batch --output verifications/salmon_merged_outs.md5 data/*/star_salmon/*.tsv
