import sys, os, glob
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
mat62 = matlist.blosum62

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
    #print "###", pdb
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
    for n in xrange(1,10):
        f = "pdbs/"+pdb+"/"+pdb+"_nohet_"+str(n)+".pdb"
        for rec in SeqIO.parse(f+".fasta", "fasta"):
            id = rec.id
            desc = rec.description
            seq = rec.seq
            break
        if ref_pdbA is None:
            aln = pairwise2.align.globalds(seq, seqA, mat62, -5, -1)[0]
            #print "DB:", f, aln[2]
            if aln[2] >= len(seq)*4.0:
                #print "SeqA found!"
                ref_pdbA = f
                ref_seqA = seq
                for i in range(len(aln[0])):
                    if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                        print pdb, "seqA", i, aln[0][i], aln[1][i]
        if ref_pdbB is None:
            aln = pairwise2.align.globalds(seq, seqB, mat62, -5, -1)[0]
            #print "DB:", f, aln[2]
            if aln[2] >= len(seq)*4.0:
                #print "SeqB found!"
                ref_pdbB = f
                ref_seqB = seq
                for i in range(len(aln[0])):
                    if aln[0][i]!=aln[1][i] and aln[0][i]!="-" and aln[1][i]!="-":
                        print pdb, "seqB", i, aln[0][i], aln[1][i]
        if ref_pdbA != None and ref_pdbB != None:
            print ref_pdbA, ref_pdbB
            break
    #check
    if ref_pdbA is None or ref_pdbB is None:
        print "A/B chain ERROR!"
        break

    #match target chain
    seqT = es[t["antigen_seq"]].strip()
    ref_seqT = ""
    ref_pdbT = None
    for n in xrange(1,10):
        f = "pdbs/"+pdb+"/"+pdb+"_nohet_"+str(n)+".pdb"
        for rec in SeqIO.parse(f+".fasta", "fasta"):
            id = rec.id
            desc = rec.description
            seq = rec.seq
            break
        if ref_pdbT is None:
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
        print "target chain ERROR!"
        break
