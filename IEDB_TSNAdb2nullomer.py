#!/usr/bin/python
# encoding: utf-8
#Produces a tab seperated file with nullomers for neoepitope in each row in the neoepitope file.
#nullomer_list  mutated_neoepitope_coding_sequence  WT-neoepitope_pair  transcript_ID  functional_annotation
#peptide IDs generated for each transcript as follows: first,last,mutated peptide wrt. original inpuut peptide
#Usage:
#python3 IEDB_TSNAdb2nullomer.py  <neoepitope files IEDB ENST2 | TSNAdb_frequent_ICGC_per_ENST3> <Homo_sapiens cds 16mers> <Homo_sapiens GRCh38 cds fasta> <output_prefix>
#Example:
#python3 IEDB_TSNAdb2nullomer.py TSNAdb_frequent_ICGC_per_ENST3.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_ICGC

'''
@author:     Mayank Mahajan

@copyright:  2022 Brigham and women's hospital. All rights reserved.

@license:    GPL3

@contact:    yashumayank@gmail.com
@deffield    updated: Updated

'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import numpy
import sys, csv, re
import itertools

d = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    '_': ['TAA', 'TAG', 'TGA'],
}

def generator(protein):
    l = [d[aa] for aa in protein]
    for comb in itertools.product(*l):
        yield comb

DELETION, INSERTION, MATCH = range(3)
def smith_waterman(seq1, seq2, insertion_penalty = -1, deletion_penalty = -1,
                   mismatch_penalty = -1, match_score = 2):
    """
    Find the optimum local sequence alignment for the sequences `seq1`
    and `seq2` using the Smith-Waterman algorithm.
    """
    m, n = len(seq1), len(seq2)

    # Construct the similarity matrix in p[i][j], and remember how
    # we constructed it -- insertion, deletion or (mis)match -- in
    # q[i][j]

    p = numpy.zeros((m + 1, n + 1))
    q = numpy.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            deletion = (p[i - 1][j] + deletion_penalty, DELETION)
            insertion = (p[i][j - 1] + insertion_penalty, INSERTION)
            if seq1[i - 1] == seq2[j - 1]:
                match = (p[i - 1][j - 1] + match_score, MATCH)
            else:
                match = (p[i - 1][j - 1] + mismatch_penalty, MATCH)
            p[i][j], q[i][j] = max((0, 0), deletion, insertion, match)

    # Yield the aligned sequences one character at a time in reverse
    # order.
    def backtrack():
        i, j = m, n
        while (i > 0 or j > 0):
            assert i >= 0 and j >= 0
            if q[i][j] == MATCH:
                i -= 1
                j -= 1
                yield seq1[i], seq2[j]
#                yield seq1[i], seq2[j], str(int(op-p[i][j]))
#                op=p[i][j]
            elif q[i][j] == INSERTION:
                j -= 1
                yield '-', seq2[j]
#                yield '-', seq2[j], str(int(op-p[i][j]))
#                op=p[i][j]
            elif q[i][j] == DELETION:
                i -= 1
                yield seq1[i], '-'
#                yield seq1[i], '-', str(int(op-p[i][j]))
#                op=p[i][j]
            else:
                assert(False)
    return [s[::-1] for s in zip(*backtrack())]

kmerPep=set()
wtEpi={}
pepLen={}
neoEpi={}
failA, failD, failI, failS, failS0 = 0, 0, 0, 0, 0
null_len = 0
padding = 6
flank = 15

print("Reading epitope db")
#read all the neo-epitopes together with (WT/frameshift) ENST, plen DB_name
#concatenate DB_names (order of reading: TCGA, ICGC, IEDB, IEDB_tcellbcell, dbNeoPep)
with open(sys.argv[1]) as IEDBf:
    for lines in IEDBf:
        IEDBtok = lines.split()
        if(IEDBtok[0] in wtEpi):
            wtEpi[IEDBtok[0]] = wtEpi[IEDBtok[0]] + ";" + IEDBtok[3]
            neoEpi[IEDBtok[0]] = neoEpi[IEDBtok[0]] + ";" + IEDBtok[1]
            pepLen[IEDBtok[0]] = pepLen[IEDBtok[0]] + ";" + IEDBtok[2]
            print("inside wtEpi -if-: detected duplicate ENST")
        else:
            wtEpi[IEDBtok[0]] = IEDBtok[3]
            neoEpi[IEDBtok[0]] = IEDBtok[1]
            pepLen[IEDBtok[0]] = IEDBtok[2]

print("Reading cds k-mers")
with open(sys.argv[2]) as kmerCdna:
    for kmerSeq in kmerCdna:
        kmerTok = kmerSeq.split()
        kmerPep.add(kmerTok[0])
        null_len = len(kmerTok[0])

nf = open(sys.argv[4] + "_neoepitopes-nullomers.tsv", "w")

#--------------- generic script to account for upto 2 SNVs per neoepitope. ----------------
print("Searching cds for wt-epitopes and derive nullomers from the respective neoepitopes")
for rec in SeqIO.parse(sys.argv[3], "fasta"):
    recIDtok = rec.id.split(".")
    try:
        ann = rec.description.split("scription:")[1].split(" [")[0]
    except IndexError:
#        print(rec.description)
        ann = "Unannotated in Ensembl"
    if(recIDtok[0] in wtEpi):
#        print(recIDtok[0] + "\t" + wtEpi[recIDtok[0]])
        for frame in range(3):
            aa_seq=(rec.seq[frame:].translate(to_stop=False))
            wtEpiTok= wtEpi[recIDtok[0]].split(";")
            neoEpiTok= neoEpi[recIDtok[0]].split(";")
            for wti in range(len(wtEpiTok)):
                try:
                    wtStart = (aa_seq.index(wtEpiTok[wti]) * 3)  + frame
                except ValueError:
                    continue

                wtEnd = (len(wtEpiTok[wti]) * 3) + wtStart
#                print(str(wtEpiTok[wti])  +"->"+ str(neoEpiTok[wti]) + ":")
#                print(rec.seq[wtStart:wtEnd].translate(to_stop=False))
                try:
                    wt,neo =smith_waterman("Z"*padding + wtEpiTok[wti] + "Z"*padding, "Z"*padding + neoEpiTok[wti] + "Z"*padding)
                except AssertionError:
                    print(str(wtEpiTok[wti])  + "->" + str(neoEpiTok[wti]) + ":")
                    failA=failA+1
                    continue
                snpPos, snpaa = [], []
                sc = 0
                for wtidx in range(padding,len(wt)-padding):
                    if(wt[wtidx]=='-'):
#                        print("i " + str(wtidx))
#                        print(wtEpiTok[wti]  +"->"+ neoEpiTok[wti])
                        failI=failI+1
                        sc=3
                        break
                    elif(neo[wtidx]=='-'):
#                        print("d " + str(wtidx))
#                        print(wtEpiTok[wti]  +"->"+ neoEpiTok[wti])
                        failD=failD+1
                        sc=3
                        break
                    elif(wt[wtidx]!=neo[wtidx]):
#                        print(neo[wtidx] + " s " + str(wtidx-padding))
                        snpPos.append(wtStart + ((wtidx-padding) * 3))
                        snpaa.append(neo[wtidx])
                        sc=sc+1
                        if(sc>2):
                            failS=failS+1
                            print(wtEpiTok[wti]  +"->"+ str(neoEpiTok[wti]))
                            break
                if(sc>2):
                    continue

                sp1=wtStart - flank
                st1=wtEnd + flank
                if(sp1<0):
                    sp1=0
                if(st1>len(rec.seq)):
                    st1=len(rec.seq)
                if(sc==1):
                    neoCds_list = [rec.seq[sp1:snpPos[0]] + lc + rec.seq[snpPos[0]+3:st1] for lc in (d[snpaa[0]])]
                elif(sc==2):
                    print(wtEpiTok[wti]  +"->"+ str(neoEpiTok[wti]) + "\t" + str(snpPos[0])+";"+str(snpPos[1]))
                    neoCds_list = [rec.seq[sp1:snpPos[0]] + lc1 + rec.seq[snpPos[0]+3:snpPos[1]] + lc2 + rec.seq[snpPos[1]+3:st1] for lc1,lc2 in generator("".join(snpaa))]
                else:
                    print(wtEpiTok[wti]  +"->"+ str(neoEpiTok[wti]) + "\t" + ';'.join(str(snpPos)))
                    failS0=failS0+1
                    continue
                for neocds in neoCds_list:
                    laststt=len(neocds)-null_len
                    null_str = ""
#                    stt_str, stp_str = "", ""
                    for stt in range(laststt):
                        stp = stt + null_len
                        if not(neocds[stt:stp] in kmerPep):
                            null_str = str(neocds[stt:stp]) + ";" + null_str
#                            stt_str = str(wtStart+stt+1) + ";" + stt_str
#                            stp_str = str(wtStart+stp) + ";" + stp_str
                    if(len(null_str) != 0):
                        nf.write(null_str[:-1] + "\t"+ str(neocds) +"\t"+ str(wtEpiTok[wti])  +"->"+ str(neoEpiTok[wti]) + "\t" + str(recIDtok[0]) +  "\t" + ann + "\n")

print("A:" + str(failA) + " D:" + str(failD) +" I:" + str(failI) +" S:" + str(failS) +" S0:" + str(failS0))
nf.close()
