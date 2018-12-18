import matplotlib
import re
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import style, colors
from matplotlib import gridspec

def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def exclusive(a, b):
    """ return items either exclusive in a or b """
    return list(set(a) ^ set(b))

def unique2a(a, b):
    """ return items only in a, but not in b """
    return list(set(a) & set(set(a) ^ set(b)))

def list2csv(lst, of):
    with open(of, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        for val in lst:
            writer.writerow([val]) 

def parse_strelka_vcf(vcf):
    df = pd.read_csv(vcf, comment='#', sep='\t', header=None, low_memory=False)
    patient = vcf.split('/')[4]
    df = df[[0,1,3,4,7]]
    df.columns = ['chr', 'pos', 'ref', 'alt', 'effect']
    df = df[(df['effect'].str.contains("HIGH"))|(df['effect'].str.contains("MODERATE"))|(df['effect'].str.contains("LOW"))]
    if not df.empty:
        df['impact'], df['impact_type'], df['gene'] = df['effect'].apply(lambda x: parse_effect(x)).str.split('@').str
        df['patient'] = patient
    df = df.drop('effect', axis=1)
    return df

# keep it easy for now pick HIGH, MODERATE and then LOW
def parse_effect(line):
    effs = line.split('EFF=')[1].split(',')
    #     extract impact, impact_type and gene
    effs = ['@'.join(list(np.array(re.split('\(|\|',ef))[[0,1,6]])) for ef in effs if ('HIGH' in  ef) or ('MODERATE' in ef) or ('LOW' in ef)]
    effs = list(set(effs))
    high = [ef for ef in effs if 'HIGH' in ef]
    moderate = [ef for ef in effs if 'MODERATE' in ef]
    low = [ef for ef in effs if 'LOW' in ef]
    if high:
        anno = high[0]
    elif moderate:
        anno = moderate[0]
    elif low:
        anno = low[0]
    else:
        print('ERROR!')
#     make sure the genes have the same name
    genes = [ef.split('@')[2] for ef in effs]
    
    return anno

