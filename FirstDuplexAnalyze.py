import os
#os.chdir('C:\\Users\\user\\Documents\\thesis\\Tools\\GiladIsana')
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pickle
from Bio import SeqIO
from collections import Counter
import RNA
from RichDuplex import *
from SeedFeatures import *



import seaborn as sns


human_clash_data_utr3 = pd.read_csv("Data/Human/Parsed/human_clash_data_utr3.csv")

biomart_df = pd.read_csv("Data/Human/Parsed/3utr.csv")

for i in range (100) :

    mirna = human_clash_data_utr3.miRNA_seq[i]
    mrna = human_clash_data_utr3.mRNA_seq_extended[i]
    if human_clash_data_utr3.ensg[i]=="ENSG00000137309":
        print "debug"
    dp = RichDuplex (mirna, mrna)
    print "ENSG: {} ENST: {} MIR: {}".format(human_clash_data_utr3.ensg[i],human_clash_data_utr3.enst[i],human_clash_data_utr3.microRNA_name[i])
    print (dp)
    c_seed = dp.IRP.extract_seed()
    print c_seed
    seed_mrna, seed_mir = dp.IRP.seed_without_bulges(c_seed)
    compact_seed = InteractionRichPresentation ("", seed_mrna, seed_mir,"")
    seed_feature = SeedFeatures (compact_seed)
    seed_feature.extract_seed_features()
    print "CLASH file seed type: {}".format(human_clash_data_utr3.seed_type[i])
    print seed_feature



    print "mrna_seed =  {} ".format(seed_mrna)
    print "mirna_seed = {} ".format(seed_mir)

    print ("**********************************************************************")
print ("gilad")