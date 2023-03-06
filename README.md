# rsApi

1. Description

Retrieve SNP annotation data from NCBI's dbSNP API from a given list of rs numbers.

2. Requirements

Python modules:

-requests

-re

-sys

3. Usage

rsApi.py listOfRSnums

The given "listOfRSnums" should be a text file with one RefSNP number per line. E.g.:

rs110
rs2293930
rs188440
(...)

