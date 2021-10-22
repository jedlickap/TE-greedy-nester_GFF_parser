#! /usr/bin/python3

import sys
import re
import json
from collections import namedtuple

gff = sys.argv[1]

# gff ="gff_for_test/Papt_JHAMT_nester_out.gff"
# gff_file to dict
# dict structure
# {'te_id':{'seq':'seq_name','te_base':'int', 'feature':{'nested_repeat': ['start,'end','attributes'], 'repeat_fragment': [[]], 
# 'polypeptide_conserved_region': [[],[]], 'primer_binding_site': [], 'RR_tract': [], 'long_terminal_repeat': [[], []]}}}}
def feature_list(l):
    l = l.rstrip()
    index_list = [3,4,8]    # indexes for start, end and attributes in gff line
    fList = [l.split("\t")[i] for i in index_list]
    return fList

def fill_in_dict(gff):
    attr_list = ["repeat_fragment", "polypeptide_conserved_region", 
    "primer_binding_site", "RR_tract", 
    "long_terminal_repeat", "target_site_duplication"]
    in_dict = {}
    te_id = 0
    reg = r'=(TE_BASE \d+);'
    chrom = ""
    te = ""
    with open(gff) as f:        
        for l in f:
            if "##" not in l:                    
                if "nested_repeat" in l:
                    teID = "TE_" + str(te_id)
                    chrom = str(l.split("\t")[0])
                    te = re.search(reg, l).group(1).replace("TE_BASE ", "")
                    in_dict[teID] = {'seq':chrom, 'te_base':te,'feature':{'nested_repeat':[]}}
                    in_dict[teID]['feature']["nested_repeat"] = feature_list(l)
                    te_id += 1
                elif chrom in l and te + ";" in l and "nested_repeat" not in l:
                    for attr in attr_list:
                        if attr in l:
                            if attr not in in_dict[teID]['feature'].keys():
                                in_dict[teID]['feature'][attr] = []
                                in_dict[teID]['feature'][attr].append(feature_list(l))
                            else:
                                in_dict[teID]['feature'][attr].append(feature_list(l))                
    return in_dict

# dict2object transformation:
# source: https://www.geeksforgeeks.org/python-convert-string-to-tuple/
class obj:
    # constructor
    def __init__(self, dict1):
        self.__dict__.update(dict1)

def dict2obj(dict1):
    # using json.loads method and passing json.dumps
    # method and custom object hook as arguments
    return json.loads(json.dumps(dict1), object_hook=obj)       
# iterate through gff obj class:
# [r for r in dir(te_obj.TE_0.feature) if not r.startswith('__')]
# analysis steps:
# gff to dict
# dict to obj
# work with obj via 'dir()' iterations and using filters (e.g. domains) and/or subsequent analysis e.g. LTR similarity
# iterate through object attributes and supply string back to attributes:
# for t in [r for r in dir(te_obj) if not r.startswith('__')]:
#     te_obj_t = (getattr(te_obj, t))
#     print(te_obj_t.seq)

# gff line parser:
# gff_line example: ['2457', '3228', 'ID=DOMAIN 0-0-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM GAG_crm;name=GAG;color=#225ea8']
class gff_line:
    
    # constructor
    def __init__(self, list1):
        self.start = int(list1[0])
        self.end = int(list1[1])
        self.attribute = list1[2]
        atr_list = [i.split("=")[1] for i in list1[2].split(";")]
        if "Parent=" not in list1[2]:
            self.ID = atr_list[0]
            self.name = atr_list[1]
            self.color = atr_list[3]
        elif "annot=" in list1[2]:
            self.ID = atr_list[0]
            self.Parent = atr_list[1]
            self.annot = atr_list[2]
            self.name = atr_list[3]
            self.color = atr_list[4]
        else:
            self.ID = atr_list[0]
            self.Parent = atr_list[1]
            self.name = atr_list[2]
            self.color = atr_list[3]

"""
TODO: solve fragmented TEs
TODECIDE: using OOP or pure functional and/or procedural notation
Next steps:
A. Modul for counting basic stuff (including generation of TE IDs lists):
    1. seq count; TEs count - fragmented versus non-fragmented
    2. filter all TEs with TSD -> tsdTEs
    3. number of tsdTEs annotated 
        3.1 at least one prot. domain annotated
    4. tsdTEs -> filt for presence of PBS and PPT -> ppTsdtTEs
        4.1 ppTsdtTEs with prot. domain (at least one)
        4.2 ppTsdtTEs with NO prot. domain -> potential LARDs and/or
B. Modul for generation of GFFs based on filtered TE IDs lists
    
"""
