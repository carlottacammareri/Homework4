#exercise1grph

from Bio import SeqIO                                                  
prefixes = []                                                          
suffixes = []                                                          
handle = open('rosalind_grph.txt', 'r')                                 
for record in SeqIO.parse(handle, 'fasta'):                            
    count1 = 0                                                         
    count2 = 0                                                         
    prefix = [record.id]                                               
    suffix = [record.id]                                               
    pre = ''                                                           
    suf = ''                                                           
    for nt in record.seq:                                              
        if count1 < 3:                                                 
            pre += nt                                                  
            count1 += 1                                                
    prefix.append(pre)                                                 
    for tn in reversed(record.seq):                                    
        if count2 < 3:                                                 
            suf += tn                                                  
            count2 += 1                                                
    suffix.append(''.join(reversed(suf)))                              
    prefixes.append(prefix)                                            
    suffixes.append(suffix)                                            
handle.close()                                                         
                                                                       
for i, k in enumerate(suffixes):                                       
    currentsf = suffixes[i][1]                                         
    currentid = suffixes[i][0]                                         
    for j, l in enumerate(prefixes):                                   
        if currentsf == prefixes[j][1] and currentid != prefixes[j][0]:
            print(currentid, prefixes[j][0])


#exerise2mprt
from bs4 import BeautifulSoup
import requests
import os
import re

def get_html(url):
   _html = ""
   resp = requests.get(url)
   if resp.status_code == 200:
      _html = resp.text
   return _html

id, seq, start, total =[], [], [], []

f = open("rosalind_mprt.txt","r")
lines = f.readlines()
for line in lines:
   id.append(line.replace("\n", ""))
f.close()
p = re.compile(r"N[^P](S|T)[^P]")
for ide in id :
   URL = "https://www.uniprot.org/uniprot/%s.fasta" % ide
   html = get_html(URL)
   soup = BeautifulSoup(html, 'html.parser')
   f = open('rosalind_mprt.txt' + "temp.txt", "w")
   f.write(soup.getText())
   f.close()
   f = open('rosalind_mprt.txt' + "temp.txt", "r")
   lines = f.readlines()
   seq= ""
   totemp=[]
   for j in range(1,len(lines)):
       seq = seq + lines[j].replace("\n","")

   for k in range(len(seq)+1):
        a = seq[k:]
        m = p.match(a)
        if m is not None:
           totemp.append(k+1)
   if totemp == [[]]:
      total.append('None')
   else:
      total.append(totemp)
   f.close()
   total = list(map(str,total))
for z in range(len(id)):
   if total[z] == '[]':
      continue
   else:
      print(id[z])
      print(total[z].replace("[","").replace("]","").replace(",",""))
os.remove('rosalind_mprt.txt'+"temp.txt")

#exercise3orf

dna_codon_table = {
    "TTT":"F", "CTT":"L", "ATT":"I", "GTT":"V",
    "TTC":"F", "CTC":"L", "ATC":"I", "GTC":"V",
    "TTA":"L", "CTA":"L", "ATA":"I", "GTA":"V",
    "TTG":"L", "CTG":"L", "ATG":"M", "GTG":"V",
    "TCT":"S", "CCT":"P", "ACT":"T", "GCT":"A",
    "TCC":"S", "CCC":"P", "ACC":"T", "GCC":"A",
    "TCA":"S", "CCA":"P", "ACA":"T", "GCA":"A",
    "TCG":"S", "CCG":"P", "ACG":"T", "GCG":"A",
    "TAT":"Y", "CAT":"H", "AAT":"N", "GAT":"D",
    "TAC":"Y", "CAC":"H", "AAC":"N", "GAC":"D",
    "TAA":"STOP", "CAA":"Q", "AAA":"K", "GAA":"E",
    "TAG":"STOP", "CAG":"Q", "AAG":"K", "GAG":"E",
    "TGT":"C", "CGT":"R", "AGT":"S", "GGT":"G",
    "TGC":"C", "CGC":"R", "AGC":"S", "GGC":"G",
    "TGA":"STOP", "CGA":"R", "AGA":"R", "GGA":"G",
    "TGG":"W", "CGG":"R", "AGG":"R", "GGG":"G"
}
def complementary_dna(dna_string):
    replace_bases = {"A":"T","T":"A","G":"C","C":"G"}
    return ''.join([replace_bases[base] for base in reversed(dna_string)])
def translate(fragment):
    peptide = []
    methionine = fragment.find("ATG")
    codons = [fragment[methionine:methionine+3] for methionine in range(methionine, len(fragment), 3)]
    for codon in codons:
        peptide += dna_codon_table.get(codon)
    return ''.join(map(str, peptide)) 
import re
 
pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAA|TAG|TGA))')
 
fragments = []
dna_string = 'TTAACGTACATTTTAACTTTTTCTCAAGCATCGTGCTTGGGGGTGCGGACCCGCCGGTTTACAATCCCTACTGATGCAGCATATTAGGAGGGTCCCATTCTCATGGCATAGGAGATTGGTCGCCCTTACGGTTCTCCAAGATCTCCTGTAATTAATCGAGCTCTTAGAACAATTTAATGTTGCGAGGCAGAAAGTATATACAACTGTATCGTGGGGGACGATATATGGACCTAGCGGCAAATGGCCCTAGGAACCTGCAGCTGTAGTGAATTCACACAGTATAGTCTGTGCCCCTCTCTCTCCCATTTTCAGCCAGTTATAAGCACCGTCCCTGCTTTTTGATCCCCCTGGCTGAGTGGCGCGATAAGCCCGCTCACAATGCCGTCTGCCACCCACGCCGAAAGTTTCGCTTATGTATAATCTATCTCATGTGATGTTTTAAGGCACGGGAACGTAGTAGCTACTACGTTCCCGTGCCTTAAAACATGCGAGGAACCTTCAGCCCATTATAGAGCTGTTATACAATTAGTCTAGGAAGGACCATTCGGAGACATGATGTACATCTACACGGACACAGAAACAGCACTCATTTCAGGATCGAACTAAGATGATCTCTGCCCTGACAGCACAGATTGTTAGCACATATAACATTGGCCCATGTTCCCTCCAAATCGCAGCGAACATAAGGCCTCCTCACTACGCTAGGTAAACTACATTCGGCCAACATATTAAAGACACGCTTGATACTTGCTAAGATGATAGCTGCCAGATCTTTGCTTTCTTTATCTAACATCTGACAGTTGTCTGGGTCGGAGTAAGTCGAAATAAGGAACTACCCCCTGAAAATCAGACAAAGACCAGAATCTGTAGCTACGTATCTTAGTAACGTAGAGACGGTTAATATGGTACGGGATGCACACG'
for string in re.findall(pattern, dna_string):
    fragments.append(string)
 
for string in re.findall(pattern, complementary_dna(dna_string)):
    fragments.append(string)
with open("rosalind_orf.txt") as file:
    for line in file:
        if line.startswith(">"):
            nextline = str()
        else:
            nextline += (line.strip("\n"))
    dna_string = nextline
 
for string in set(fragments):
    print(translate(string))

#exercise4splc

RNA_CODON_TABLE = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#read input file and generate the mRNA
def processInput():
    inputFile = open('rosalind_splc.txt','r')
    sequences = inputFile.readlines()
    i=1
    tempString=''
    inputStrings =[]
    while(i<len(sequences)):
        if(sequences[i].__contains__('>')):
            inputStrings.append(tempString)
            tempString=''
            i+=1
        elif(i==len(sequences)-1):
            inputStrings.append(sequences[i])
            i+=1
        else:
            tempString +=sequences[i]
            i+=1
    inputStrings[0] = inputStrings[0].replace('\n', '')
    i=1
    while (i < len(inputStrings)):
        inputStrings[i] = inputStrings[i].replace('\n','')
        inputStrings[0] = removeIntrons(inputStrings[0],inputStrings[i])
        i+=1
    inputStrings[0] = inputStrings[0].replace('T','U')
    return inputStrings[0]


def removeIntrons(dnaSequence, intron):
    dnaSequence = dnaSequence.replace(intron, '')
    return dnaSequence


def replaceCodons(rnasequence):
    i = 0
    protein = ""
    rnasequence = rnasequence.replace("\n","")
    while (i < len(rnasequence)):
        if(RNA_CODON_TABLE[rnasequence[i:i+3]]== 'STOP'):
            break
        protein += RNA_CODON_TABLE[rnasequence[i:i + 3]]
        i += 3
    protein = protein.replace('STOP','')
    print(protein)

replaceCodons(processInput())


#exercise5pmch

from Bio import SeqIO                      
from math import factorial                 
sequence = ''                              
handle = open('rosalind_pmch.txt', 'r')     
for record in SeqIO.parse(handle, 'fasta'):
    sequence = str(record.seq)             
handle.close()                             

AU = 0                                     
GC = 0                                     
for nt in sequence:                        
    if nt == 'A':                          
        AU += 1                            
    elif nt == 'G':                        
        GC += 1                            

matchings = factorial(AU) * factorial(GC)  
print(matchings)   


#exercise6pper

n = 80
k = 8
partial_perm = 1
for i in range(k):
    partial_perm *= (n - i)
print(partial_perm % 1000000)

#exercise7tree

data = []                                          
with open('rosalind_tree.txt', 'r') as f:             
    for line in f:                                 
        split_data = [int(x) for x in line.split()]
        data.append(split_data)                    

n = data[0][0]                                     
edges = data[1:]                                   
print(n - len(edges) - 1) 

#exercise8long

f = open("rosalind_long.txt", "r")
mat = []
str1 = f.read()
str1 = str1.replace("Rosalind_", "")
str1 = str1.replace("\n", "")
str1 = ''.join([i for i in str1 if not i.isdigit()])
mat = str1.split(">")
mat.remove("")

def long(mat, bbb=''):
    if (len(mat) == 0):
        return bbb

    elif (len(bbb) == 0):
        bbb = mat.pop(0)
        return long(mat, bbb)

    else:
        for i in range(len(mat)):
            a = mat[i]
            for j in range(len(a) // 2):
                c = len(a) - j
                if bbb.startswith(a[j:]):
                    mat.pop(i)
                    return long(mat, a[:j] + bbb)
                if bbb.endswith(a[:c]):
                    mat.pop(i)
                    return long(mat, bbb + a[c:])
print(long(mat))

#exercise9sseq

def find_motif(data, motif):
    position, indices = -1, ''
    for nucleotide in motif:
        position = data.find(nucleotide, position+1)
        indices += str(position+1) + ' '
    print(indices)


with open('rosalind_sseq.txt', 'r') as file:
    content = file.read()
DNAs_number, lines, line_number, DNAs = content.count('>'), content.splitlines(), 0, []
for i in range(DNAs_number):
    DNA = ''
    line_number += 1
    while lines[line_number][0] != '>':
        DNA += lines[line_number]
        line_number += 1
        if line_number+1 > len(lines):
            break
    DNAs.append(DNA)

find_motif(DNAs[0], DNAs[1])

#exercise10tran

from Bio import SeqIO

f = open("rosalind_tran.txt", 'r')
raw_data = list(SeqIO.parse(f, "fasta"))
f.close()

seq1 = raw_data[0].seq
seq2 = raw_data[1].seq

transitions = 0
transversions = float(0)
for nt in range(len(seq1)):
    if seq1[nt] == seq2[nt]:
        continue
    elif seq1[nt] == "A" and (seq2[nt] == "C" or seq2[nt] == "T"):
        transversions += 1
    elif seq1[nt] == "G" and (seq2[nt] == "C" or seq2[nt] == "T"):
        transversions += 1
    elif seq1[nt] == "C" and (seq2[nt] == "A" or seq2[nt] == "G"):
        transversions += 1
    elif seq1[nt] == "T" and (seq2[nt] == "A" or seq2[nt] == "G"):
        transversions += 1
    else:
        transitions += 1

print (transitions/transversions)


