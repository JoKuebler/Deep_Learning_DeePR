#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 08:57:38 2017

Creates a JSON representation of a HHR file. 
Also allows pretty printing of that HHR file such that it can be
better printed into a file

@author: Lukas Zimmermann
"""
from __future__ import division,print_function
import sys
import re
import json
from collections import defaultdict
from itertools import chain


#%% TOKEN DEFINITIONS
JSONSEP = ','
OBJOPEN = '{'
OBJCLOSE = '}'
QUOTE = '"'
ARROPEN = "["
ARRCLOSE = "]"            

#%% dict to JSON 

def dict2JSONString(d):
    return "".join([OBJOPEN,
                    JSONSEP.join(('"{}": {}'.format(key, py2JSONString(value)) 
                    for (key,value) in d.items())),
                    OBJCLOSE])
# %% List to JSON


def list2JSONString(d):
    return "".join([ARROPEN,
                    JSONSEP.join((py2JSONString(value) for value in d)),
                    ARRCLOSE])
#%% Something else to JSON    
    
# Represents a Python object as JSON String, if possible
def py2JSONString(d):
    if isinstance(d, int):
        return str(d)
    if isinstance(d, float):
        return str(d)
    if isinstance(d, str):
        return '{}{}{}'.format(QUOTE, d, QUOTE)
    if isinstance(d, dict) or isinstance(d, defaultdict):
        return dict2JSONString(d)
    if isinstance(d, list):
        return list2JSONString(d)
#%%% Regular Expressions used for parsing the HHR file
info_line = re.compile(r"\s*(\S+)\s+(\S.*)$")
hit_line = re.compile(r"""
([0-9]+)\s+
(\S+)\s+
(.*)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
([0-9Ee\.\+\-]+)\s+
([0-9Ee\.\+\-]+)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
([0-9]+)\s+
([0-9]+)-([0-9]+)\s+
([0-9]+)-([0-9]+)\s*
\(([0-9]+)\)\s*
$
""".replace('\n',''))

# Second version of hitline, where struc is not explicitly available
hit_line2 = re.compile(r"""
([0-9]+)\s+
(.*)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
([0-9Ee\.\+\-]+)\s+
([0-9Ee\.\+\-]+)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
(-?[0-9]+(?:\.[0-9]+)?)\s+
([0-9]+)\s+
([0-9]+)-([0-9]+)\s+
([0-9]+)-([0-9]+)\s*
\(([0-9]+)\)\s*
$
""".replace('\n',''))


# Terminates the alignment
align_term = re.compile(r"\s*[Nn]o\s+([0-9]+)\s*")


#%% Returns the metainfo section of the HHR file as a PyDict
    
# Returns a dictionary of the Info Part of the HHR file
# With the input iterator
def fetch_info(it):
    res = {}
    for line in it:
        line = line.strip()
        # Stop at first empty line
        if not line:
            break
        grps = re.match(info_line, line).groups()
        res[grps[0]] = grps[1]
    return res
    



#%% Returns the Hitlist of the HHR file as an array of dicts

# Annotates header and corresponding python type
headers = [('no',int), 
          ('struc', str),
          ('hit', str),
          ('prob', float),
          ('eval', float),
          ('pval', float),
          ('score', float),
          ('ss', float),
          ('cols', int),
          ('query_begin', int),
          ('query_end', int),
          ('template_begin', int),
          ('template_end', int),
          ('ref', int)]

# Alternative version of headers, where struc is not identifiable
headers2 = [('no',int), 
            ('struc', str),
            ('prob', float),
            ('eval', float),
            ('pval', float),
            ('score', float),
            ('ss', float),
            ('cols', int),
            ('query_begin', int),
            ('query_end', int),
            ('template_begin', int),
            ('template_end', int),
            ('ref', int)]

def fetch_hits(it):
    res = []
    state = 0
    for line in it:
        line = line.strip()
        if not line and state == 1:
            break
        if not line and state == 0:
            continue
        # Trash the header line
        if line.startswith('No') or line.startswith('no'):
            state = 1
            continue
        # Try to match hit lines
        m = re.match(hit_line, line)
        h = headers
        if m is None:
           m = re.match(hit_line2, line)
           h = headers2
           if m is None:
               print("Error matching hitline:\n" + line)
               sys.exit()
        res.append({ head[0] : head[1](grp)  
        for (head, grp) in  zip(h, m.groups())})
    return res
        
        
    
#%% Returns the Alignments of the HHR file as an array of dicts
#VALVA; sum_probs was originally float; however in the HHR file
#is is often Inf. Since we don't parse this value for use by the
#Toolkit, str should handle all cases

header_align_info = [ ("probab", float),
                      ("eval", float),
                      ("score", float),
                      ("aligned_cols", int),
                      ("identities", lambda x: float(x.replace("%",''))/100.0),
                      ("similarity", float),
                      ("sum_probs", str),
                      ("template_neff", float)]


# Probab=99.72  E-value=1.6e-18  Score=173.01  Aligned_cols=197  Identities=13%  Similarity=0.094  Sum_probs=0.0

# Fetches the info line of one alignment and returns it as hash
def readinfo(line):
    return {head[0]: head[1](grp) for (head, grp) in zip(header_align_info,
                                       (pair.split("=")[1] 
                                       for pair in line.split())) }
    

    

ks = ["ss_pred","ss_dssp", "consensus", "seq"]


# Probab=99.66  E-value=5e-17  Score=159.50  Aligned_cols=195  Identities=14%  Similarity=0.102  Sum_probs=0.0


crange = None



def handlesequenceline(line, d):
    global crange
    spt = line.split()
    word = spt[0]
    # Handle Secondary Structure Lines
    if word == "ss_pred":
        d["ss_pred"]+=spt[1]
    elif word == "ss_dssp":
        d["ss_dssp"] += spt[1]
    elif word == "Consensus":
        d["consensus"]+=spt[2]
        if d["start"] == -1:
            d["start"] = int(spt[1])
        d["end"] = int(spt[-2])
        d["ref"] = int(spt[-1][1:-1])
    # We must be in the line containing the name
    else:
        tmpName = spt[0]
        if len(tmpName.split("|")) > 1:
             tmpName = tmpName.split("|")[1]
        d["name"] = tmpName
        d["seq"] += spt[2]
        q = line.index(spt[2])
        crange = (q+2, q+2+len(spt[2]))
        if d["start"] == -1:
            d["start"] = int(spt[1])
        d["end"] = int(spt[-2])
        d["ref"] = int(spt[-1][1:-1])
    return d
            

# Fetches one particular alignment of the alignments list as JSON
def parse_alignment(it):
    query = defaultdict(list)
    template = defaultdict(list)
    res = defaultdict(list)
    query["consensus"] = ""
    query["end"] = -1
    query["name"] = ""
    query["ref"] = -1
    query["seq"] = ""
    query["ss_pred"] = ""
    query["ss_dssp"] = ""
    query["start"] = -1
    template["consensus"] = ""
    template["end"] = -1
    template["struc"] = ""
    template["ref"] = -1
    template["seq"] = ""
    template["ss_dssp"] = ""
    template["ss_pred"] = ""
    template["start"] = -1
    stop_line = None

    number_seen = False
    for line in it:
        line = line.rstrip('\n')
        # Skip empty lines
        if not line.strip():
            continue
        # Set number of alignment or determine whether the alignment is complete
        if not number_seen and (line.startswith("No") or line.startswith("no")):
            number_seen = True
            res["no"] = re.match(align_term, line).groups()[0]         
            continue
        if number_seen and (line.startswith("No") or line.startswith("no")):
            stop_line = line
            break
        
        if line.startswith("Confidence"):
            res["confidence"].append(line[crange[0]:crange[1]])
            
        elif line.startswith('>'):
            res["header"] = ' '.join(line.split()[1:])
        elif line.startswith("P") or line.startswith("p"):
            res["info"] = readinfo(line)
        elif line.startswith("Q") or line.startswith("q"):
            query = handlesequenceline(re.sub(r"^Q\s+","",line), query)
        elif line.startswith("T"):
            template = handlesequenceline(re.sub(r"^T\s+","",line), template)
        else:
            res["agree"].append(line[crange[0]:crange[1]])       
    for k in ks:
        if k in template:
            template[k] = ''.join(template[k])
        if k in query:
            query[k] = ''.join(query[k])
    res["query"] = query
    res["template"] = template
    res["agree"] = "".join(res["agree"])
    
    if "confidence" in res:
        res["confidence"] = "".join(res["confidence"])
    
    return stop_line, res
    
def main(argv):
    
    # HHR File represented as dictionary
    hhr = {}
    with  open(argv[1], 'r') as infile:
        hhr["info"] = fetch_info(infile)
        hhr["hits"] = fetch_hits(infile)
        stop_line = ""
        alignments = []
        it = infile
        while  stop_line is not None:
            stop_line, ali = parse_alignment(it)
            alignments.append(ali)
            # Stop line has to be reparsed
            it = chain([stop_line], it)
        hhr["alignments"] = alignments
        print(json.dumps(hhr, indent=4, sort_keys=True))

if __name__ == '__main__':
    main(sys.argv)
    
