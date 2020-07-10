#!/usr/bin/python

import argparse
import sys
import os.path
import subprocess
import itertools
import shutil
import sqlite3
import csv
def storeindexfile(indexfile):
    strain2spacer={}
#    for line in indexfile:
#       print 1
#    with open(indexfile, 'r') as spacer2strain: #strain to spacer file in PDI
    f=open('strain2spacer.csv','w')
    f.write('dataset_id,bact_strain,strain_id,spacer\n')
    x=0
    for strain in indexfile:
           strain2spacer[strain.split('\t')[0]]=strain.strip().split('\t')[1].split(' ')
           x=x+1
           for spacer in strain.strip().split('\t')[1].split(' '):
                f.write(','.join(['1',strain.split('\t')[0],'strain_'+str(x),spacer])+'\n')
    f.close()       
    return strain2spacer

def complementpam(PAMlist):
    newlist=[]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for dna in PAMlist:
        cpam=''.join([complement[base] for base in dna])
        newlist.append(cpam)
    return newlist

def makeblastdb(fasta):
    subprocess.call(['makeblastdb','-parse_seqids','-dbtype','nucl','-in', fasta])

def PAMProto(spacerfile,dbfile):
    subprocess.call(['python2', 'PAMProtoPatternGrab_full.py', spacerfile ,dbfile])

def viruslist(virusfasta):
    with open(virusfasta) as fasta:
        virL=[line.strip()[1:] for line in fasta if '>' in line]
    return virL 

def storealnfile(aln):
    phage2spacer={}
    f=open('virus2spacer.csv', 'w')
    f.write('dataset_id,virus_strain,virus_cluster,spacer,criteria\n')
    with open(aln, 'r') as alnfile:
        for line in alnfile:
             phage2spacer.setdefault(line.split('\t')[1],[]).append(line.split('\t')[0])
             f.write(','.join(['1',line.split('\t')[1],'NA',line.split('\t')[0],'NA'])+'\n')
    f.close()
    return phage2spacer

def storealnfilereads(aln):
    phage2spacer={}
    #f.open('virus2spacer.csv')
    #f.write('dataset_id,virus_strain,virus_cluster,spacer,criteria')
    with open(aln, 'r') as alnfile:
        for line in alnfile:
            phage2spacer.setdefault(aln,[]).append(line.split('\t')[0])
          #','.join('1',line.split('\t')[1],line.split('\t')[0],'NA')
    return phage2spacer
#def storealnfile(aln):
#    phage2spacer={}
#    with open(aln, 'r') as alnfile:
#        for line in alnfile:
#            phage2spacer.setdefault(line.split('\t')[1],[]).append(line.split('\t')[0])
#    return phage2spacer

def PAMfilter(PAMs, placutoff, alnfile, output):
    f=open(output,'w')
    with open(alnfile, 'r') as filteredextra:
        for line in filteredextra:
	    array=line.split('\t')
	    PLA= array[17] 
	    PAM=array[18] #pseudo PAM
		#oprint PAM
            if len(PAMs) > 0:
	      if PAM in PAMs or PAM[1:] in PAMs:
	        if float(PLA) >= float(placutoff):
                    f.write(line)
            else:
	        if float(PLA) >= float(placutoff):
                    f.write(line)
		    #print line.strip()
    f.close()

def makesqlitedb(strain2spacer, virus2spacer,databasename):
    conn=sqlite3.connect(databasename)
    cur=conn.cursor()
    cur.execute("CREATE TABLE bacteria (dataset_id,bact_strain, strain_id, spacer);")
    cur.execute("CREATE TABLE data_sets (id,name, Description, original_virus_filename, original_bact_filename,provider);")
    cur.execute("CREATE TABLE virus (dataset_id, virus_strain,virus_cluster,spacer,criteria);")
    with open(strain2spacer, 'rb') as s2s:
        dr=csv.DictReader(s2s)
        #print dr
        tdb=[(i['dataset_id'],i['bact_strain'],i['strain_id'],i['spacer']) for i in dr]
        #print tdb
    with open(virus2spacer, 'rb') as v2s:
        dr1=csv.DictReader(v2s)
        #print dr1
        tdb1=[(i['dataset_id'],i['virus_strain'],i['virus_cluster'],i['spacer'],i['criteria']) for i in dr1]
        #print tdb1
    cur.executemany("INSERT INTO bacteria (dataset_id,bact_strain, strain_id, spacer) VALUES (?,?,?,?);",tdb)
    cur.executemany("INSERT INTO virus (dataset_id,virus_strain, virus_cluster, spacer,criteria) VALUES (?,?,?,?,?);",tdb1)
    conn.commit()
    conn.close()



def runPDI(phage2spacer, strain2spacer, virL,output):
    f=open(output,'w')
    f.write('phage\tPI\tPDI\tIDI\n')
    print virL
    for phage in virL:
      totalimmune=0
      PDIcount=0
      PDI=0
      totalspacer=[]
      if phage in phage2spacer:
        Ndict={}
        spacername2count={}
        for strain in strain2spacer:
            spacername=str(' '.join(strain2spacer[strain]))
            spacer2proto=list(set(strain2spacer[strain]).intersection(set(phage2spacer[phage])))
            if len(spacer2proto)>0:
                totalimmune=totalimmune+1
                totalspacer=totalspacer+spacer2proto
            if spacername in spacername2count:
                spacername2count[spacername]=spacername2count[spacername]+1
            else:
                spacername2count[spacername]=1
        PI=float(totalimmune)/len(strain2spacer)
        IDI=len(totalspacer)/float(len(strain2spacer))
        ################relative abundance##########
        for spacername in spacername2count:
            Ndict[spacername]=spacername2count[spacername]/float(len(strain2spacer))
        comparisonlist=list(itertools.permutations(list(set(Ndict.keys())),2))
        ##################PDI calculation############
        for combo in comparisonlist:
            sigma=0
            maxN=max(Ndict[combo[1]],Ndict[combo[0]])
            PDItemp=1-abs(Ndict[combo[1]]-Ndict[combo[0]]/maxN)
            spacer2proto1=list(set(combo[0].split(' ')).intersection(set(phage2spacer[phage])))
            spacer2proto2=list(set(combo[1].split(' ')).intersection(set(phage2spacer[phage])))
            unsharedmatch1= set(spacer2proto1).difference(spacer2proto2)
            unsharedmatch2= set(spacer2proto2).difference(spacer2proto1)
            if len(unsharedmatch1) > 0 and len(unsharedmatch2) > 0:
                 sigma=Ndict[combo[0]]*Ndict[combo[1]]*(1/Ndict[combo[0]])
                 PDI=PDI+sigma*PDItemp
        #############################################
        s= '\t'.join([phage,str(PI), str(PDI),str(IDI)])+'\n'
        f.write(s)
      elif phage not in phage2spacer:
        s= '\t'.join([phage, '0','0','0'])+'\n'
        f.write(s)
    f.close()  

def runPDIread(phage2spacer, strain2spacer, output):
    f=open(output,'w')
    f.write('phage\tPI\tPDI\tIDI\n')
    #print virL
    for phage in phage2spacer:
      totalimmune=0
      PDIcount=0
      totalspacer=[]
      if phage in phage2spacer:
        for strain in strain2spacer:
            spacer2proto=list(set(strain2spacer[strain]).intersection(set(phage2spacer[phage])))
            if len(spacer2proto)>0:
                totalimmune=totalimmune+1
                totalspacer=totalspacer+spacer2proto
        PI=float(totalimmune)/len(strain2spacer)
        IDI=len(totalspacer)/float(len(strain2spacer))
        comparisonlist=list(itertools.combinations(strain2spacer,2))
        for combo in comparisonlist:
            spacer2proto1=list(set(strain2spacer[combo[0]]).intersection(set(phage2spacer[phage])))
            spacer2proto2=list(set(strain2spacer[combo[1]]).intersection(set(phage2spacer[phage])))
            unsharedmatch1= set(spacer2proto1).difference(spacer2proto2)
            unsharedmatch2= set(spacer2proto2).difference(spacer2proto1)
            if len(unsharedmatch1) > 0 and len(unsharedmatch2) > 0:
                 PDIcount=PDIcount+1
        PDI=float(PDIcount)/len(comparisonlist)
        s= '\t'.join([phage,str(PI), str(PDI),str(IDI)])+'\n'
        f.write(s)
      elif phage not in phage2spacer:
        s= '\t'.join([phage, '0','0','0'])+'\n'
        f.write(s)
    f.close()  
    
def main():
    parser = argparse.ArgumentParser(description="Do something.")
    parser.add_argument('-p', '--PAM', required=True, nargs='*')
    parser.add_argument('-c', '--cutoff', required=True)
    parser.add_argument('-s', '--spacerfasta', required=True)
    parser.add_argument('-i', '--indexfile', type=argparse.FileType('r'), default=sys.stdin,required=True)
    parser.add_argument('-v', '--virusfasta', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-q', '--complement',action='store_true')
    parser.add_argument('-r', '--reads',action='store_true')
    parser.add_argument('-d', '--sqldatabase', required=True)
    args = parser.parse_args()
    strain2spacer=storeindexfile(args.indexfile) 
    print 'PAM list:',args.PAM
    if args.complement == True:
        PAMs=complementpam(args.PAM)
        print 'complement PAM list:',PAMs
    else:
        PAMs=args.PAM
    alnfilename=os.path.basename(args.virusfasta)[:]+'_vs_'+os.path.basename(args.spacerfasta)
    extra=alnfilename+'.dir/'+alnfilename+'.extra.aln'
    directory=alnfilename+'.dir/'
    if os.path.exists(extra):
        filt=extra+'.PAMcutoff'
    else:
        print "building blast db with:", args.virusfasta
        makeblastdb(args.virusfasta)
        print "blasting and extending:", args.spacerfasta,args.virusfasta
        if os.path.isdir(directory):
            print "diretory exists, removing directory and rerunning"
            shutil.rmtree(directory)
        PAMProto(args.spacerfasta, args.virusfasta) 
        filt=extra+'.PAMcutoff'
    print "filtering blast"
    PAMfilter(PAMs, args.cutoff, extra, filt)
    if args.reads == False: 
        phage2spacer=storealnfile(filt)
        virL=viruslist(args.virusfasta)
        runPDI(phage2spacer,strain2spacer, virL,args.output)
        print "you are running complete viral genomes"
    elif args.reads == True:
        phage2spacer=storealnfilereads(filt)
        runPDIread(phage2spacer,strain2spacer, args.output)
        print "you are running reads"
    makesqlitedb('strain2spacer.csv','virus2spacer.csv',args.sqldatabase)
        
if __name__== '__main__':
    main()
