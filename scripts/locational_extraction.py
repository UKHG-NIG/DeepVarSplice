#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 10:22:20 2022

@author: emre
"""

import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table


maf_file = pd.read_csv('shared_vcf_tumor.maf',delimiter='\t')

# Define the input GTF file
gtf_file = "Homo_sapiens.GRCh38.105.gtf"

# Define a dictionary to store gene symbols and their exon coordinates
exon_dict = {}

def gene_mut(gene,fin_dat,maf):
    df=fin_dat[fin_dat['Gene']==gene]
    df.drop_duplicates(subset=['exon_start', 'exon_end'], inplace=True)
    df['length'] = df['exon_end'] - df['exon_start'] + 1
    
    #df = df[(df['exon_start']>39463000) & (df['exon_start']<39485000)]
    loc = ','.join(map(str,df['exon_start']))
    loca = gene + '_locs <- c(' + loc + ')'
    length = ','.join(map(str,df['length']))
    lengtha = gene + '_lengths <- c(' + length + ')'
    mut = ','.join(map(str,maf[maf['Hugo_Symbol']==gene]['vcf_pos']))
    muta = gene + '_muts <- c(' + mut + ')'
    print('#######' + gene + '#######')
    print(loca)
    print(lengtha)
    print(muta)


# Open the GTF file
with open(gtf_file) as f:
    # Iterate over each line in the file
    for line in f:
        # Skip comments and blank lines
        if line.startswith("#") or line.strip() == "":
            continue
        
        # Split the line into fields
        fields = line.strip().split("\t")
        
        # Skip lines that are not exons
        if fields[2] != "exon":
            continue
        
        # Extract the chromosome, start, end coordinates, and gene symbol
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        gene_id = None
        gene_name = None
        
        # Extract the gene symbol and ID from the attributes field
        attributes = fields[8].split(";")
        for attr in attributes:
            if attr == '':
                continue
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1]
            elif attr.startswith(" gene_name"):
                gene_name = attr.split('"')[1]
            elif attr.startswith(" exon_id"):
                exon_id = attr.split('"')[1]
        # Skip exons without a gene symbol
        if gene_name is None:
            continue
        
        # Add the exon coordinates to the dictionary for the corresponding gene
        if gene_name in exon_dict:
            exon_dict[gene_name].append((chrom, start, end, exon_id))
        else:
            exon_dict[gene_name] = [(chrom, start, end, exon_id)]

# Print the number of genes and the number of unique exon coordinates
print("Number of genes:", len(exon_dict))
print("Number of unique exon coordinates:", sum(len(v) for v in exon_dict.values()))

new_dict = {k: list(set(v)) for k, v in exon_dict.items()}

genes = []
chromosome = []
start = []
end = []
exon_id = []

for k in new_dict.keys():
    for v in new_dict[k]:
        genes.append(k)
        chromosome.append(v[0])
        start.append(v[1])
        end.append(v[2])
        exon_id.append(v[3])
    
final_data = pd.DataFrame({'Gene': genes,'Chromosome': chromosome, 'exon_start': start, 'exon_end':end, 'Exon_ID':exon_id})

#Gene list to use if there are specific genes to be extracted

#genes = ['RBFOX1','CSMD1','WWOX','BAGE2','SMYD3','MTUS2','CDH13','AC013287.1',
#         'MAGI2','LRP1B','DLG2','PTPRN2','EXOC4','JMJD1C','LINC02055',
#         'NRG1','FHIT','SLC30A10','SNX29','HHAT','NRXN3','NAV2','ALK',
#         'MACROD2','CNTN5','PDE1C','PRKN','NOS1AP','FRMD4A','CNTNAP2',
#         'UBBP4','THSD4','KMT2C']


gene_list=[]
total_number = []
num_ss = []
num_mid = []
num_deep = []
exon_list = []
for gene in list(set(final_data['Gene'])):
    red_maf = maf_file[maf_file['Hugo_Symbol']==gene].reset_index(drop=True)
    total = len(red_maf['Hugo_Symbol'])
    red_coords = final_data[final_data['Gene']==gene].reset_index(drop=True)
    
    red_coords['ss_1'] = red_coords['exon_start']-2
    red_coords['ss_2'] = red_coords['exon_end']+2
    
    red_coords['mid_1'] = red_coords['ss_1']-30
    red_coords['mid_2'] = red_coords['ss_2']+30
    
    total_number.append(total)
    sss = 0
    mids= 0
    exons = 0
    gene_list.append(gene)
    for i in range(len(red_coords['Exon_ID'])):
        
        strt = red_coords['exon_start'][i]
        endd = red_coords['exon_end'][i]
        
        ss_1 = red_coords['ss_1'][i]
        ss_2 = red_coords['ss_2'][i]
        
        mid_1 = red_coords['mid_1'][i]
        mid_2 = red_coords['mid_2'][i]
        
        sss1 = red_maf['Start_Position'].between(ss_1,strt).sum()
        sss2 = red_maf['Start_Position'].between(endd,ss_2).sum()
        
        sss = sss1 + sss2 + sss
        
        mid1 = red_maf['Start_Position'].between(mid_1,ss_1).sum()
        mid2 = red_maf['Start_Position'].between(ss_2,mid_2).sum()
        
        exon_1 = red_maf['Start_Position'].between(strt,endd).sum()
         
        exons = exons + exon_1
        
        mids = mid1 + mid2 + mids
        
    num_ss.append(sss)
    num_mid.append(mids)
    exon_list.append(exons)
    num_deep.append(total-(mids + sss))
        
final_result = pd.DataFrame({'Gene_ID':gene_list,
                            'Number_of_Total':total_number,
                            'Number_of_SS':num_ss,
                            'Number_of_Prox':num_mid,
                            'Number_of_Deep':num_deep})        

final_result.to_csv('Grouped_mutations_all.csv',sep='\t')


by_total = final_result.sort_values(by='Number_of_Total',ascending=False)[0:60].reset_index(drop=True)
by_ss = final_result.sort_values(by='Number_of_SS',ascending=False)[0:60].reset_index(drop=True)
by_prox = final_result.sort_values(by='Number_of_Prox',ascending=False)[0:60].reset_index(drop=True)

ax = plt.subplot(111, frame_on=False) # no visible frame
ax.xaxis.set_visible(False)  # hide the x axis
ax.yaxis.set_visible(False)  # hide the y axis
table(ax, by_prox, loc='center')  # where df is your data frame

plt.savefig('by_prox_fin.pdf',format='pdf', dpi=800,bbox_inches='tight')
 
genes = ['MACF1']


for gene in genes:
    gene_mut(gene,final_data,maf_file)


file = open('proximal_mutations.txt','w')
file = open('ss_mutations.txt','w')
file = open('deep_mutations.txt','w')

for i in list(by_prox['Gene_ID']):
    item = "'" + i + "',"
    file.write(item)
file.close()

