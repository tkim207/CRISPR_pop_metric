#!/usr/bin/python

import argparse
import sys
import os.path
import subprocess
import itertools
import shutil

def storeindexfile(indexfile):
    strain2spacer={}
#    for line in indexfile:
#       print 1
#    with open(indexfile, 'r') as spacer2strain: #strain to spacer file in PDI
    for strain in indexfile:
           strain2spacer[strain.split('\t')[0]]=strain.strip().split('\t')[1].split(' ')
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
    with open(aln, 'r') as alnfile:
        for line in alnfile:
            phage2spacer.setdefault(line.split('\t')[1],[]).append(line.split('\t')[0])
    return phage2spacer

def storealnfilereads(aln):
    phage2spacer={}
    with open(aln, 'r') as alnfile:
        for line in alnfile:
            phage2spacer.setdefault(aln,[]).append(line.split('\t')[0])
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

def runPDI(phage2spacer, strain2spacer, virL,output):
    f=open(output,'w')
    f.write('phage\tPI\tPDI\tIDI\n')
    print virL
    for phage in virL:
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
##        IDI=len(totalspacer)/float(totalimmune)
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
##        IDI=len(totalspacer)/float(totalimmune)
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
        

 
    
#    opts = ap.parse_args([])
#    if not any([opts.a, opts.b, opts.c]):
#        ap.print_usage()
#        quit()

#    print("This won't run.")

if __name__== '__main__':
    main()
