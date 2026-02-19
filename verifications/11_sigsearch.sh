#!/bin/bash

cap md5 --slurm batch --output verifications/sigsearch.md5 --ignore "*_safety_data.csv" "*.rds" ./results/signature_reversion/*
