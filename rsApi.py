#!/usr/bin/python3

__author__ = 'DFernandes'


import requests
import re
import sys

outList = []
outList.append("dbSNP_ID\tVariant_Type\tGene\tGene_Full\tPathology_Name\tPathological_Significance\tConsequence\tAminoacid_Change\tChromosome\tPosition_hg19\tPosition_hg38\tReference_Allele\tDerived_Allele\tFrequency_dbGaP_Ref\tFrequency_dbGaP_Der\tsamtools_pileup_hg19\tgatk_pileup_hg19\tsamtools_pileup_hg38\tgatk_pileup_hg38")
rsList = open(sys.argv[1], 'r').readlines()

for it in rsList:
    if(len(it) > 1):
        oldRS = ""
        rsID = it.strip()
        idNoRs = str(rsID.split('rs')[1])
        print(rsID)
        link = str("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"+idNoRs)
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36'}
        f = requests.get(link, headers=headers)

        ## If merged SNP, get data from new one
        if "merged_snapshot_data" in f.json().keys():
            oldRS = rsID
            oldIdNoRS = idNoRs
            idNoRs = f.json()['merged_snapshot_data']['merged_into'][0]
            rsID = "rs"+str(idNoRs)

            link = str("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + idNoRs)

            f = requests.get(link, headers=headers)
        ## Get variant type and reference allele
        #print(f.json()['primary_snapshot_data']['allele_annotations'][0].keys()) # find frequencies here 'frequency'
        varType = str(f.json()['primary_snapshot_data']['variant_type'])
        try:
            refAl = f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][0]['allele']['spdi']['deleted_sequence']
        except:
            refAl = " "
        ## Get all described alternative alleles and frequencies
        freqsList = {}
        derAl = []
        derAlOr = []
        ## Get alleles
        try:
            if len(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles']) > 2:
                for alu in range(1,len(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'])):
                    derAl.append(re.split('(\d+)', f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][alu]['hgvs'].split(":g.")[1])[2].split(">")[1])
                derAlOr = derAl
                derAl = "/".join(derAl)
            elif len(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles']) == 2:
                derAl = re.split('(\d+)',f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][1]['hgvs'].split(":g.")[1])[2].split(">")[1]
                derAlOr = derAl
        except:
            derAl = ""
        freqsList = {i : "" for i in derAlOr}
        ## Get global frequencies from dbGaP or 1000Genomes
        try:
            it = 0
            ale = ""
            aleFr = ""
            for i in f.json()['primary_snapshot_data']['allele_annotations']:
                try:
                    if len(list(filter(lambda person: person['study_name'] == 'dbGaP_PopFreq', (i['frequency'])))) != 0:
                        ale = list(filter(lambda person: person['study_name'] == 'dbGaP_PopFreq', (i['frequency'])))[0]['observation']['inserted_sequence']
                        aleFr = int(list(filter(lambda person: person['study_name'] == 'dbGaP_PopFreq', (i['frequency'])))[0]['allele_count']) / int(list(filter(lambda person: person['study_name'] == 'dbGaP_PopFreq', (i['frequency'])))[0]['total_count'])
                    if len(list(filter(lambda person: person['study_name'] == 'dbGaP_PopFreq', (i['frequency'])))) == 0 and len(list(filter(lambda person: person['study_name'] == '1000Genomes', (i['frequency'])))) != 0:
                        ale = list(filter(lambda person: person['study_name'] == '1000Genomes', (i['frequency'])))[0]['observation']['inserted_sequence']
                        aleFr = int(list(filter(lambda person: person['study_name'] == '1000Genomes', (i['frequency'])))[0]['allele_count']) / int(list(filter(lambda person: person['study_name'] == '1000Genomes', (i['frequency'])))[0]['total_count'])
                    if ale != "":
                        freqsList[ale] = aleFr
                except:
                    pass
                it += 1
            majorFreq = round(freqsList[refAl],5)
        except:
            majorFreq = ""
            minorFreq = ""
            minorFreqAl = ""
        ## Transform alleles and frequencies into text format for output
        freqsList2 = []
        derAl = []
        for i1 in freqsList:
            if i1 != refAl and len(str(freqsList[i1])) != 0:
                freqsList2.append(str(i1)+"="+str(round(freqsList[i1],5)))
            if i1 != refAl:
                derAl.append(i1)
        freqsList2 = "/".join(freqsList2)
        derAl = "/".join(derAl)
        if len(str(majorFreq)) > 0:
            majorFreq = refAl + "=" + str(majorFreq)
        ## If SNV and no frequencies exist, get alternative allele from somewhere else
        if varType == "snv" and len(refAl) != 0 and len(derAl) == 0:
            derAl = (f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][1]['hgvs'].split(":g.")[1].split(">")[1])  # 158630656A>G
        ## Get gene name, clinical description, and significance
        try:
            geneAcro = str(f.json()['primary_snapshot_data']['allele_annotations'][0]['assembly_annotation'][0]['genes'][0]['locus'])
            geneLong = str(f.json()['primary_snapshot_data']['allele_annotations'][0]['assembly_annotation'][0]['genes'][0]['name'])
        except:
            geneAcro = ""
            geneLong = ""
        try:
            clinName = str(f.json()['primary_snapshot_data']['allele_annotations'][1]['clinical'][0]['disease_names'][0]) # Last [0] might restrict, in case there are more than one name
            clinSign = str(f.json()['primary_snapshot_data']['allele_annotations'][1]['clinical'][0]['clinical_significances'][0]) # Last [0] might restrict, in case there are more than one name
            aminoChange = str(f.json()['primary_snapshot_data']['placements_with_allele'][int(len(f.json()['primary_snapshot_data']['placements_with_allele']))-1]['alleles'][1]['hgvs'].split(':p.')[1])
        except:
            clinName = ""
            clinSign = ""
            aminoChange = ""
        if len(clinName) == 0:
            try:
                clinName = str(f.json()['primary_snapshot_data']['allele_annotations'][2]['clinical'][0]['disease_names'][0])
            except:
                clinName = ""
        try:
            conseq = str(f.json()['primary_snapshot_data']['allele_annotations'][1]['assembly_annotation'][0]['genes'][0]['rnas'][0]['protein']['sequence_ontology'][0]['name'])
        except:
            conseq = ""

        #print(f.json()['primary_snapshot_data']['placements_with_allele'][1]) ## [0] = hg38, [1] = hg19
        #print(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles']) # alleles and positions
        #print(len(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'])) # gives total number of alleles
        #print(f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][1]['hgvs'].split(":g.")[1]) # 158630656A>G
        ## Get genomica position in hg19 and hg38
        chromo = re.split('0+', f.json()['primary_snapshot_data']['allele_annotations'][0]['assembly_annotation'][0]['seq_id'])[1].split('.')[0]
        if chromo == '12920':
            chromo = "MT"
        hg19pos = re.split('(\d+)',f.json()['primary_snapshot_data']['placements_with_allele'][1]['alleles'][1]['hgvs'].split(":g.")[1])[1]
        hg38pos = re.split('(\d+)',f.json()['primary_snapshot_data']['placements_with_allele'][0]['alleles'][1]['hgvs'].split(":g.")[1])[1]
        ## Append to list of SNPs to output at end
        if len(oldRS) == 0:
            outList.append(rsID+"\t"+varType+"\t"+geneAcro+"\t"+geneLong+"\t"+clinName+"\t"+clinSign+"\t"+conseq+"\t"+aminoChange+"\t"+chromo+"\t"+hg19pos+"\t"+hg38pos+"\t"+refAl+"\t"+derAl+"\t"+str(majorFreq)+"\t"+str(freqsList2)+"\t"+chromo+" "+str(int(hg19pos)-1)+" "+str(hg19pos)+"\t"+chromo+":"+str(hg19pos)+"-"+str(hg19pos)   +"\t"+chromo+" "+str(int(hg38pos)-1)+" "+str(hg38pos)+"\t"+chromo+":"+str(hg38pos)+"-"+str(hg38pos))
        elif len(oldRS) > 0:
            outList.append(oldRS + "_mergedInto_" + rsID + "\t" + varType + "\t" + geneAcro + "\t" + geneLong + "\t" + clinName + "\t" + clinSign + "\t" + conseq + "\t" + aminoChange + "\t" + chromo + "\t" + hg19pos + "\t" + hg38pos + "\t" + refAl + "\t" + derAl + "\t" + str(majorFreq) + "\t" + str(freqsList2) + "\t" + chromo + " " + str(int(hg19pos) - 1) + " " + str(hg19pos) + "\t" + chromo + ":" + str(hg19pos) + "-" + str(hg19pos) + "\t" + chromo + " " + str(int(hg38pos) - 1) + " " + str(hg38pos) + "\t" + chromo + ":" + str(hg38pos) + "-" + str(hg38pos))

    #print(rsID+"\t"+varType+"\t"+"\t"+chromo+"\t"+hg19pos+"\t"+hg38pos+"\t"+refAl+"\t"+derAl+"\t"+str(refAl)+"="+str(majorFreq)+"\t"+str(freqsList2)+"\t"+chromo+" "+str(int(hg19pos)-1)+" "+str(hg19pos)+"\t"+chromo+":"+str(hg19pos)+"-"+str(hg19pos)   +"\t"+chromo+" "+str(int(hg38pos)-1)+" "+str(hg38pos)+"\t"+chromo+":"+str(hg38pos)+"-"+str(hg38pos))
fout = open(sys.argv[1] + "_rsAPIoutput.txt",'w')
for ite in outList:
    fout.write(ite)
    fout.write("\n")
fout.close()
