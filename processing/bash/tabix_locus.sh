#! /usr/bin/bash

## 1) Input tsv file
## 2) Output file name

bgzip -f $1;
tabix -f -S1 -b 2 -e 2 -s 3 $2
