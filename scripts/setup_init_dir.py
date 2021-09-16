import sys, os, glob
import os.path
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from collections import defaultdict
mat62 = matlist.blosum62
contact_cut = 8

#map light/heavy chain
templates = []
for rec in SeqIO.parse("template.fasta", "fasta"):
    templates.append( rec.seq )
assert( len(templates) == 2 )

def check_A_B_order( A, B ):
    aln = pairwise2.align.globalds(templates[0], A, mat62, -5, -1)[0]
    scL = aln[2]
    aln = pairwise2.align.globalds(templates[1], B, mat62, -5, -1)[0]
    scH = aln[2]
    return scL>500.0 and scH>500.0

lines = open(sys.argv[1]).readlines()
t = {}
for n, tag in enumerate(lines[0].strip().split()):
    t[tag] = n
#print t

for l in lines[1:]:
    es = l.strip().split('\t')
    pdb = es[0]
    #check interface
    map_partner = defaultdict(list)
    map_partner_small = defaultdict(list)
    fi = "pdbs/" + pdb + "/interface.txt"
    for ll in open(fi, 'r').readlines():
        tmp = ll.strip().split()
        c = int(tmp[2])
        a = int(tmp[0])
        b = int(tmp[1])
        map_partner_small[a].append(b)
        map_partner_small[b].append(a)
        if c >= contact_cut:
            map_partner[a].append(b)
            map_partner[b].append(a)
    #print map_partner

    seqA = es[t["antibody_seq_a"]].strip()
    seqB = es[t["antibody_seq_b"]].strip()
    #switch
    #if seqA[:3] in ["EVQ", "QVQ", "QVK", "GVQ", "VQL", "EIS", "SEV"] and seqB[:3] in ["DIV", "EIV", "DIQ", "DIE", "VLT", "AIV", "EIA", "ALT"]:
    #    seqA, seqB = seqB, seqA
    if not check_A_B_order( seqA, seqB ):
        seqA, seqB = seqB, seqA
    dG = float(es[4])

    #match AB chain
    ref_seqA = ""
    ref_seqB = ""
    ref_pdbA = None
    ref_pdbB = None
    ref_idA = 0
    ref_idB = 0
    for n in xrange(1,10):
        f = "pdbs/"+pdb+"/"+pdb+"_nohet_"+str(n)+".pdb"
        #check exist
        if not os.path.isfile(f): break
        for rec in SeqIO.parse(f+".fasta", "fasta"):
            id = rec.id
            desc = rec.description
            seq = rec.seq
            break
        if ref_pdbA is None:
            flag = True
            if ref_pdbB != None:
                if n not in map_partner[ref_idB]:
                    flag = False
            if flag:
                aln = pairwise2.align.globalds(seq, seqA, mat62, -5, -1)[0]
                #print "DB:", f, aln[2]
                if aln[2] >= len(seq)*4.0:
                    #print "SeqA found!"
                    ref_pdbA = f
                    ref_seqA = seq
                    ref_idA = n
                    for i in range(len(aln[0])):
                        if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                            print pdb, "seqA", i, aln[0][i], aln[1][i]
        if ref_pdbB is None:
            flag = True
            if ref_pdbA != None:
                if n not in map_partner[ref_idA]:
                    flag = False
            if flag:
                aln = pairwise2.align.globalds(seq, seqB, mat62, -5, -1)[0]
                #print "DB:", f, aln[2]
                if aln[2] >= len(seq)*4.0:
                    #print "SeqB found!"
                    ref_pdbB = f
                    ref_seqB = seq
                    ref_idB = n
                    for i in range(len(aln[0])):
                        if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                            print pdb, "seqB", i, aln[0][i], aln[1][i]
        if ref_pdbA != None and ref_pdbB != None:
            break
    #check
    if ref_pdbA is None and ref_pdbB is None:
        print "ERROR: A/B chain ERROR!"
        continue

    #match target chain
    seqT = es[t["antigen_seq"]].strip()
    ref_seqT = ""
    ref_pdbT = None
    for n in xrange(1,10):
        f = "pdbs/"+pdb+"/"+pdb+"_nohet_"+str(n)+".pdb"
        if not os.path.isfile(f): break
        for rec in SeqIO.parse(f+".fasta", "fasta"):
            id = rec.id
            desc = rec.description
            seq = rec.seq
            break
        if ref_pdbT is None:
            flag = False
            if ref_pdbA != None:
                if n in map_partner[ref_idA]:
                    flag = True
            if ref_pdbB != None:
                if n in map_partner[ref_idB]:
                    flag = True
            if flag:
                aln = pairwise2.align.globalds(seq, seqT, mat62, -5, -1)[0]
                #print "DB:", f, aln[2]
                if aln[2] >= len(seq)*2.0:
                    #print "SeqT found!"
                    ref_pdbT = f
                    ref_seqT = seq
                    #print format_alignment(*aln)
                    for i in range(len(aln[0])):
                        if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                            print pdb, "seqT", i, aln[0][i], aln[1][i]
                    break
    #check
    if ref_pdbT is None:
        contact_cut = 5
    else:
        print ref_pdbA, ref_pdbB, ref_pdbT
        continue

    for n in xrange(1,10):
        f = "pdbs/"+pdb+"/"+pdb+"_nohet_"+str(n)+".pdb"
        if not os.path.isfile(f): break
        for rec in SeqIO.parse(f+".fasta", "fasta"):
            id = rec.id
            desc = rec.description
            seq = rec.seq
            break
        if ref_pdbT is None:
            flag = False
            if ref_pdbA != None:
                if n in map_partner_small[ref_idA]:
                    flag = True
            if ref_pdbB != None:
                if n in map_partner_small[ref_idB]:
                    flag = True
            if flag:
                aln = pairwise2.align.globalds(seq, seqT, mat62, -5, -1)[0]
                #print "DB:", f, aln[2]
                if aln[2] >= len(seq)*2.0:
                    #print "SeqT found!"
                    ref_pdbT = f
                    ref_seqT = seq
                    #print format_alignment(*aln)
                    for i in range(len(aln[0])):
                        if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                            print pdb, "seqT", i, aln[0][i], aln[1][i]
                    break

    #check
    if ref_pdbT is None:
        print "ERROR: target chain missing!"
    else:
        print ref_pdbA, ref_pdbB, ref_pdbT

