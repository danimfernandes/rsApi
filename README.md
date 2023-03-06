# rsApi

### 1. Description

Retrieve SNP annotation data from NCBI's dbSNP API from a given list of rs numbers.

### 2. Requirements

Python modules:

-requests

-re

-sys

### 3. Usage

    rsApi.py listOfRSnums

The given "listOfRSnums" should be a text file with one RefSNP number per line. E.g.:

    rs387906590
    rs121912679
    rs387906589
    (...)

### 4. Output

A tab-spaced text file named "listofRSnums_rsAPIoutput.txt" will be created, with the following information:

    dbSNP_ID	Variant_Type	Gene	Gene_Full	Pathology_Name	Pathological_Significance	Consequence	Aminoacid_Change	Chromosome	Position_hg19	Position_hg38	Reference_Allele	Derived_Allele	Frequency_dbGaP_Ref	Frequency_dbGaP_Der	samtools_pileup_hg19	gatk_pileup_hg19	samtools_pileup_hg38	gatk_pileup_hg38
    rs387906590	snv	ACVR1	activin A receptor type 1	Progressive myositis ossificans	pathogenic	missense_variant	Arg375Pro	2	158617532	157761020	C	G			2 158617531 158617532	2:158617532-158617532	2 157761019 157761020	2:157761020-157761020
    rs121912679	snv	ACVR1	activin A receptor type 1	Progressive myositis ossificans	pathogenic	missense_variant	Gly356Asp	2	158617589	157761077	C	T			2 158617588 158617589	2:158617589-158617589	2 157761076 157761077	2:157761077-157761077
    rs387906589	snv	ACVR1	activin A receptor type 1	Brainstem glioma	likely-pathogenic	missense_variant	Gly328Val	2	158622516	157766004	C	A/T			2 158622515 158622516	2:158622516-158622516	2 157766003 157766004	2:157766004-157766004
