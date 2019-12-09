#!/usr/bin/python

import pandas as pd
import os
import re
import multiprocessing as mul
import time 
import random
import collections as cl
import itertools as it
import Queue
import sys
import getopt
import numpy as np
import gc
import subprocess

manager = mul.Manager()

lock = mul.Lock()
lock1= mul.Lock()
lock2= mul.Lock()

allfiles=os.listdir(os.getcwd())
if 'assemblies' in allfiles:	
	wkfolder='/home/q5476572/data/assemblies/'
else:
	wkfolder='/home/michalsakin/data/'

print wkfolder
savefolder='/home/q5476572/data/'
	
bowtiefolder='/home/q5476572/data/refindex/'
	
chromfolder='/home/q5476572/data/chroms/'
	
exonfolder='/home/q5476572/data/exon/'
	
genefile='/home/q5476572/data/knownGenefix.txt'
	
inputfile='all'

inputchrom='big'

unfinished=0
me=0
mc=0
ma=0
mh=0
mt=0
mb=0
md=0
ms=0
opts,args=getopt.getopt(sys.argv[1:],"ubhatdesc:f:")
for op, value in opts:
	if op=='-f':
		inputfile=value
	if op=='-c':
		inputchrom=value
	if op=='-u':
		unfinished=1
	if op=='-b':
		mb=1
	if op=='-h':
		mh=1
	if op=='-t':
		mt=1
	if op=='-a':
		ma=1
	if op=='-d':
		md=1
	if op=='-e':
		me=1
	if op=='-s':
		ms=1	


try:
	os.system('echo "" >  runlog.txt')
except:
	pass

try:
	os.system('echo "" > errorlog.txt')
except:
	pass

try:
	os.system('echo "" > memrecord.txt')
except:
	pass

def deepdel(x):
        if type(x)==type([]):
                while x!=[]:
                        if type(x[0])==type([]):
                                x.extend(x[0])
                        del x[0]
        del x

def checkstoage():

	f=open('runlog.txt',mode='a')
	for i in dir():
		try:
			f.write(str(i)+" "+str( eval(i).nbytes) +"\n")
		except:
			f.write(str(i)+' '+str( sys.getsizeof(eval(i)))+'\n' )
	f.close()

def checkmem(x=''):

	mem=os.popen("free").read()
	mem=re.findall(r'[0-9]+',mem.split('\n')[1])
	mem=int(mem[1])
	with open ('memrecord.txt',mode='a') as f:
		f.write(x+' '+str(mem)+'\n')
	f.close



def fasterfind(a,b):
	
	keya=sorted(range(len(a)), key= lambda x: a[x])
	keyb=sorted(range(len(b)), key= lambda x: b[x])
	a1=[a[x] for x in keya]
	b1=[b[x] for x in keyb]
	notin=[]
	posi=0
	result=[0 for x in a]
	for i,x in enumerate(a1):
		
		try:
			posi=b1[posi:].index(x)+posi
			result[keya[i]]=keyb[posi]
			
		except:

			result[keya[i]]=-1
			
	return result


def cleansame(r0):
	
	if r0==[]:
		return []

	r1=[r0[0]]
	for a in r0[1:]:
		
		if len(r1)>0 and a==r1[-1]:
			del r1[-1]
			continue

		r1.append(a)
		
	r0=r1
	return r0

def cutabsrange(start,end,posi):

	if posi==[]:
		return []
	
	allposi=[start,end]+posi
	
	key=sorted(range(len(allposi)), key=lambda x: allposi[x])
	
	cutindexstart=key.index(0)
	cutindexend=key.index(1)
	
		
	if cutindexend==cutindexstart+1:
		
		if cutindexstart/2*2!=cutindexstart:
			
			return [start,end]
		
		else:
			
			return []
		
	else:

		if key[cutindexstart+1]/2*2==key[cutindexstart+1]:
			
			cutindexstart=cutindexstart+1
		
		if key[cutindexend-1]/2*2!=key[cutindexend-1]:
			
			cutindexend=cutindexend-1	

		
		if allposi[key[cutindexstart]]==allposi[key[cutindexstart+1]]:
			
			cutindexstart=cutindexstart+2
		
		if allposi[key[cutindexend-1]]==allposi[key[cutindexend]]:
			
			cutindexend=cutindexend-2
			
	
		return [allposi[a] for a in key[cutindexstart:cutindexend+1]]
		


def cutrange(cutrange0,exons,exone):

	cutrange=cutrange0
	for a in xrange(len(cutrange)/2):

		cutrange[a*2+1]=cutrange[a*2+1]-1

	qposi=exone[0]-exons[0]
	posifix=exons[0]
	cutindex=0
	newposi=[]

	
	while cutindex<len(cutrange) and qposi > cutrange[cutindex]:
		newposi.append(cutrange[cutindex]+posifix+cutindex-cutindex/2*2)
		cutindex=cutindex+1

	if cutindex/2*2!=cutindex:
		
		newposi.append(exone[0])
		
	
	for i in xrange(1,len(exone)):

		posifix=posifix+exons[i]-exone[i-1]

		qposi=qposi+exone[i]-exons[i]

		if cutindex/2*2!=cutindex:

			newposi.append(exons[i])


		while cutindex<len(cutrange) and  qposi > cutrange[cutindex]:

			
			newposi.append(cutrange[cutindex]+posifix+cutindex-cutindex/2*2)
			cutindex=cutindex+1

		if cutindex/2*2!=cutindex:

			newposi.append(exone[i])


	for a in xrange(cutindex, len(cutrange)):
		
		if cutindex/2*2!=cutindex:
			newposi.append(exone[-1])

		newposi.append(cutrange[a]+posifix+cutindex-cutindex/2*2)




	return  newposi
		
	


def findoverlap(start,end):
	l=len(start)
	if l<2:
		return start,end
	allposi=start+end
	
	key0=sorted(range(2*l), key=lambda i: allposi[i])
	c=0
	s=0
	start0,end0=[],[]
	for x in key0:
		
		if x < l:
			c=c+1
			if c==1:
				s=start[x]
		else:
			c=c-1
			
		if c==0:
			start0.append(s)
			end0.append(end[x-l])
	return start0, end0


def findcommon(range1,range2):

	start=[range1[2*k] for k in xrange(len(range1)/2)]+[range2[2*k] for k in xrange(len(range2)/2)]
	
	end=[range1[2*k+1] for k in xrange(len(range1)/2)]+[range2[2*k+1] for k in xrange(len(range2)/2)]
	l1=len(range1)/2
	l=len(start)
	if l<2:
		return [],[]
	allposi=start+end
	
	key0=sorted(range(2*l), key=lambda i: allposi[i])
	c1=0
	c2=0
	s=0
	range0=[]
	for x in key0:
		
		if x < l1:
			c1=c1+1
			if c1==1 and c2>0:
				s=start[x]
		elif x>=l1 and x<l:
			c2=c2+1
			if c1>0 and c2==1:
				s=start[x]
		elif x>=l and x<l+l1:
			c1=c1-1
			if c1==0 and c2>0:
				range0.append(s)
				range0.append(end[x-l])
		else:
			c2=c2-1
			
			if c2==0 and c1>0:
				range0.append(s)
				range0.append(end[x-l])

	range00=[]
	for k in xrange(len(range0)/2):
		if 	range0[2*k]!=range0[2*k+1]:
			range00.append(range0[2*k])
			range00.append(range0[2*k+1])
	
	return range00
	



def findcodi(seq,posi):
	if posi==[]:
		return []

	if ':' in seq:
		
		reposi=re.findall(r'\d+:\d+', seq)
		
		for reposi0 in reposi:
			reposi0s=int(reposi0.split(':')[0])
			reposi0e=int(reposi0.split(':')[1])
			seq=seq.replace(reposi0,''.join(['N' for a in xrange(reposi0e-reposi0s)]))


	
	gapssize=[a.span()[1]-a.span()[0] for a in re.finditer(r'-+',seq)]
	piecesize=[a.span()[1]-a.span()[0] for a in re.finditer(r'[A-Za-z]+',seq)]
	first=0
	if seq[0]=='-':
		first=1	
	s0=0
	index=0
	gaps=[]
	posifix=[a for a in posi]
	for i in xrange(len(piecesize)):
		s0=s0+piecesize[i]
		while index<len(posi) and posi[index]<s0:
			posifix[index]=posi[index]+sum(gapssize[:i+first])
			index=index+1
	if index<len(posi):
		for i in xrange(index,len(posi)):
			posifix[index]=posi[index]+sum(gapssize)				
	return posifix

def findposi1(posi0,start=0,end=-1):
	
	mrnasize=0
	posi=[start]
	for i in xrange(len(posi0)/2):
		
		mrnasize=mrnasize+posi0[2*i+1]-posi0[2*i]
		
		if mrnasize<start or (end!=-1 and mrnasize>=end):
			continue
		
		posi.append(mrnasize)

		posi.append(mrnasize)
	
	if len(posi)==1:
		return []
	elif end==-1:
		
		return posi[:-1]
	else:
		return posi+[end]
	


def hit(i,reffile,chrom,hitstart=0):
	
	
	def findhit(sam):
		
		try:
			t=pd.read_csv(sam,sep='\t',names=range(20),header=None,dtype='str')
		except:
			return {}
		
		if len(t)<1:
			return {}
		
		t=t[t.isnull().T.any().T ==False]

		if len(t)<1:
                        return {}


		exon=[str(a).split(' ')[0] for a in list(t[0].values.flatten())]
		
		piece=list(t[2].values.flatten())
		
		posi=[int(a) for a in list(t[3].values.flatten())]
		
		output=cl.defaultdict(list)
		
		for i,exon0 in enumerate(exon):
			
			output[exon0].append(str(piece[i])+':'+str(posi[i]-1))
		output=dict(output)

		return output

	def analyzehit(record,reffile,chrom,genename,sort,genes,genee,exons,exone,strand):
		
		try:
			os.mkdir('result/{:s}'.format(chrom))
		except:
			pass
		
		output0=[]
		output1=[]
		genelist=[]
		
		start=int(genes)
		
		end=int(genee)
		
		size=start-end
		
		exons0=record.keys()
		exonp=record.values()
		
		allhits=[';'.join(a) for a in exonp]
		
		output0=[[genename,a, exons0[a],allhits[a]] for a in xrange(len(exons0))]
		
		chromfind=cl.defaultdict(list)
		
		for exon in exons0:
			for a in record[exon]:
			
				chromfind[a.split(':')[0]].append([exon,int(a.split(':')[1])])
		
		s=0
		allhits=[]
		
		for a in chromfind.keys():
			
			b=sorted(chromfind[a],key=lambda x: x[1])
				
			last1=b[0][1]
			last0=b[0][0]
			for x in b:
				if x[0]!= last0 and x[1]-last1<30 :
					allhits.append('{:s}:{:d}'.format(a,x[1]))
					s=1
				last1=x[1]
				last0=x[0]
					
		
		if s==0:
			print 'find only original gene\n'
			return []	
		
		
		if len(allhits)>0:
		
			output1.append([sort,genename,start,end,';'.join(allhits),exons,exone,strand])
			
			genelist.append(genename)
		
		f=pd.DataFrame.from_records(output0).to_csv('result/{:s}/{:s}_exonout.csv'.format(reffile,chrom),sep='\t',header=None,mode="a")
		f=pd.DataFrame.from_records(output1).to_csv('result/{:s}/{:s}_geneout.csv'.format(reffile,chrom),sep='\t', header=None,mode="a")

		return [genename]

	
	start0=time.clock()
	
	genefile=exonfolder+'{:s}_gene.csv'.format(chrom.split('.fa')[0])
	
	exonfile=exonfolder+'{:s}_exon.fasta'.format(chrom.split('.fa')[0])
	
	try:
		geneinfo=pd.read_csv(genefile,skiprows=0,sep='\t',header=None)
	except:
		return []
	
	if hitstart==0:
	
		try:
			os.system('mv result/{:s}/{:s}_exonout.csv result/{:s}/{:s}_exonout.csv_bak'.format(reffile,chrom,reffile,chrom))
		except:
			pass
		
		try:
			os.system('mv result/{:s}/{:s}_geneout.csv result/{:s}/{:s}_geneout.csv_bak'.format(reffile,chrom,reffile,chrom))
		except:
			pass


	print 'start run reference #{:d} {:s} with chrom {:s}, has genes {:d}\n'.format(i,reffile,chrom,len(geneinfo))
	
	genename=list(geneinfo[2].values.flatten())
	genes=list(geneinfo[4].values.flatten())
	genee=list(geneinfo[5].values.flatten())
	strand=list(geneinfo[3].values.flatten())
	exonsort=list(geneinfo[6].values.flatten())
	exonsort=[[int(a) for a in b.split(',')] for b in exonsort]
	exons=list(geneinfo[7].values.flatten())
	exone=list(geneinfo[8].values.flatten())
	
	with open (exonfile, mode='r') as f:
		read=f.read()
	f.close
		
	read=read.split('>')
	del read[0]
	
	size0=[len(''.join(a.split('\n')[1:])) for a in read]

	e=-1
	ss0=10000
	for j,genename0 in enumerate(genename):
		
		if j <hitstart:
			continue
		
		try:
			os.mkdir('temp/{:s}/{:s}'.format(reffile,chrom))
		except:
			pass
		
		sort0=[a for a in exonsort[j] if a>e]
		
		genes0file=savefolder+'temp/{:s}/{:s}/{:s}.fa'.format(reffile,chrom,str(j/ss0))
		
		write=''
		if sort0 != []:
		
			e=max(sort0)
			write='\n>'.join([read[a] for a in sort0 if size0[a] > 18])		
			if write != '':
				write='>'+write+'\n'
		m='a'
		if j==j/ss0*ss0:
			m='w'

		with open(genes0file, mode=m) as f2:
			
			f2.write(write)
		
		f2.close
	
	del read
	
	genelist0=[]
	exonrecord=cl.defaultdict(list)
	s=[]
	for j0,genename0 in enumerate(genename):
		
		if j0<hitstart:
			continue
		
		m=j0/ss0
		
		print '------start run reference #{:d} {:s} with chrom {:s}, run gene {:d}\n'.format(i,reffile,chrom,m)
		
		mRNAfile=savefolder+'temp/{:s}/{:s}/{:s}.fa'.format(reffile,chrom,str(m))	
		outputfile='temp/{:s}/{:s}/output_{:d}'.format(reffile,chrom,m)
		
		if j0==j0/ss0*ss0:
			
	
			start1=time.clock()
				
			cm='bowtie2 -a -f --threads 16 --mm  --very-sensitive --dovetail --mp 5,4 --ignore-quals --score-min G,0,0  --no-unal --no-hd --no-sq -x {:s} -U {:s} -S {:s}'.format('refindex/'+reffile+'_bowtie', mRNAfile ,outputfile)
			s.append(3)
			while s[-1]>0:
				try:
					os.system(cm)
					s[-1]=0
				except:
					s[-1]=s[-1]-1
					with open('errorlog.txt', mode='a') as f:
						f.write('\n repeat {:s} #{:d}'.format(reffile+'_'+str(j0),4-s[-1]))
					f.close
			
			hit0=findhit(outputfile)
			
			for exon in hit0.keys():
				
				exonrecord[exon]=hit0[exon]
			
                        end1=time.clock()
                
                        with open('runlog.txt',mode='a') as f:
                                f.write('\n {:s} {:s} {:s}---end hit, used time:{:f} \n'.format(reffile,chrom,str(m),end1-start1))
                        f.close
				
			try:
				os.system('rm {:s}'.format(mRNAfile))
			except:
				print ('unable remove {:s}\n'.format(mRNAfile))
				
			try:
				os.system('rm {:s}'.format(outputfile))
			except:
				print ('unable remove {:s}\n'.format(outputfile)) 
		
		
		
		generecord0=cl.defaultdict(list)
		
		for exon in exonsort[j0]:
			
			generecord0[str(exon)]=exonrecord[str(exon)]
		
		genelist00=analyzehit(generecord0,reffile,chrom,genename0,j0,genes[j0],genee[j0],exons[j0],exone[j0],strand[j0])
		
		genelist0=genelist0+genelist00
	
	
	end0=time.clock()
	
	print 'used time:{:f}'.format(end0-start0)
	
	return genelist0

def fixseq(reffile,title,rnaseq):

	retrofile='./result/{:s}_retro.csv'.format(reffile)
	
	try:
		t=pd.read_csv(retrofile,sep=',',skiprows=1,header=None)
	except:
		return 1
	
	seq1=list(t[len(list(t.iloc[0]))-1].values.flatten())
	chrom=[str(a) for a in list(t[12].values.flatten())]
	
	for i,seq01 in enumerate(seq1):
		print i
		
		if ':' in seq01:
			rnaseq0=rnaseq[title.index(chrom[i])]
			reposi=re.findall(r'\d+:\d+', seq01)
			
			for reposi0 in reposi:
				reposi0s=int(reposi0.split(':')[0])
				reposi0e=int(reposi0.split(':')[1])
				seq01=seq01.replace(reposi0,rnaseq0[reposi0s:reposi0e])
			
			t[len(list(t.iloc[0]))-1][[i]]=seq01+''

			del seq01
		
		
	output=t.to_csv(retrofile,sep=',',header=None,index=False,mode="w")

	return 0 



def readlastz(lastzname,reffile,c0,chrom,gene,start, end , exons ,exone,strand, piecesize,searchs,searche,lread,cycle):
	

	def findposi0(exons,exone):
		
		mrnasize=0
		posi=[0]
		for i in xrange(len(exons)):
			
			mrnasize=mrnasize+int(exone[i])-int(exons[i])
			
			posi.append(mrnasize)
		
		
		posi=posi[:-1]
		posi.sort()
		
		return mrnasize, posi
	
	
	def recoverinfor(chrom0,rstart0,rend0,size1):
		

		pieces=[]
		chrom1=[0 for a in chrom0]
		rstart1=[0 for a in rstart0]
		rend1=[0 for a in rend0]

		for i in xrange(len(chrom0)):
			
			refstart=int(re.search(r'(?<=_)\d+(?=(:|;))',chrom0[i]).group())
			pieces.append(refstart)
			rstart1[i]=rstart0[i]+refstart
			rend1[i]=rend0[i]+refstart
			chrom1[i]=chrom0[i].split('_')[0]
			
		
		return chrom1,rstart1,rend1,pieces
			

	def findcombine(start,end,qs,qe,chrom,strand,k,query,ref):

		def findshift(query0,ref0,refposi,x1,x2,y1,y2):
			if refposi<0:
				return 0
			pieceposi=findcodi(ref0,[refposi])[0]
			queryposi=pieceposi-query[:pieceposi].count('-')
			return refposi-queryposi+x1-y1-y1+y2

		
		
		l=len(start)
		if l<2:
			return [[0]],[[]]
		
		samechrom=cl.defaultdict(list)
		
		for i in xrange(l):
			
			samechrom[str(chrom[i])+strand[i]].append(i)
		
		index,gaps,combined=[],[],[]
		for c in samechrom.keys():
			
			gethers=[int(start[a]) for a in samechrom[c]]
			key0= sorted(range(len(gethers)), key=lambda x: gethers[x])		
			getheri=[[samechrom[c][a]] for a in key0]
			combined=[]
			gethere=[int(end[samechrom[c][a]]) for a in key0]
			gethers=[gethers[a] for a in key0]
			qs0=[int(qs[samechrom[c][a]]) for a in key0]
			
			gap=[[] for a in samechrom[c]]
			for j in xrange(len(samechrom[c])):
				ss=0
				for a in xrange(j+1,len(samechrom[c])):
				
					dis=max(gethere[a], gethere[j])-min(gethers[j],gethers[a])-(gethere[a]-gethers[a])-(gethere[j]-gethers[j])
					if dis<=30:
						if gethers[j]<=gethers[a]:
							bigger=a
							smaller=j
							
							
						else:
							bigger=j
							smaller=a
							

						if dis<30 and (dis > -30 or gethere[bigger]<=gethere[smaller]  or findshift(query[samechrom[c][smaller]],ref[samechrom[c][smaller]],gethere[smaller]-gethers[bigger]-1,gethers[bigger],qs[samechrom[c][bigger]],gethere[smaller],qe[samechrom[c][smaller]])==0):
							getheri[a]=getheri[j]+getheri[a]
							gethers[a]=min(gethers[j],gethers[a])
							gethere[a]=max(gethere[j],gethere[a])
							gap[a]=gap[j]+[dis]
							ss=1
							break
				if ss==0:
					index.append(getheri[j])
					gaps.append(gap[j])
		


		return index,gaps
	
	
	def cordicombine(qstart,qend,rstarto,rendo,qtext,rtext,strand):

		if len(rstarto)<2:

			return qstart,qend,rstarto,rendo,qtext,rtext


		if strand=='-':
			rstart=[-a for a in rendo]
			rend=[-a for a in rstarto]
		else:
			rstart=rstarto
			rend=rendo
		
		key1=sorted(range(len(rstart)), key=lambda x: rstart[x])
		rend=[rend[a] for a in key1]
		qstart=[qstart[a] for a in key1]
		qend=[qend[a] for a in key1]
		qtext=[qtext[a] for a in key1]
		rtext=[rtext[a] for a in key1]
		rstart=[rstart[a] for a in key1]
		
		qstart0=[qstart[0]]
		qend0=[qend[0]]
			
		e0=rend[0]
		qe0=qend[0]
		qtext0=qtext[0]
		rtext0=rtext[0]
		for i in xrange(1,len(rstart)):
			
			if rend[i]<=e0:
				continue	
			if qend[i]-qstart[i]<=10:
				continue	
			lap00=e0-rstart[i]
			if lap00<=0:
				shift=0
				newpposi=0
				newqposi=qstart[i]
			else:
				newpposi=findcodi(rtext[i],[lap00-1])[0]+1
				newqposi=newpposi-qtext[i][:newpposi].count('-')+qstart[i]
				shift=0

			


			if shift==0:
				
				gapqq=''
				gapqr=''
				gaprq=''
				gaprr=''
				if lap00<0:
	
					gaprq=''.join(['-' for a in xrange(rstart[i]-e0)])
					
					gaprr='{:d}:{:d}'.format(e0,rstart[i])


				s00=newpposi

				if newqposi-qe0>0 and newqposi-qe0<10:
				
					gapqq='{:d}:{:d}'.format(min(newqposi,qe0),max(newqposi,qe0))
					gapqr=''.join(['-' for a in xrange(newqposi-qe0)])
					qend0[-1]=qend[i]

				elif newqposi-qe0<=0 and newqposi-qe0>-10:
					gapqq=''.join(['-' for a in xrange(qe0-newqposi)])
					s00=s00+qe0-newqposi
					qend0[-1]=qend[i]
				else:
				
					qstart0.append(newqposi)
					qend0.append(qend[i])
				e0=rend[i]	
				qe0=qend[i]
				qtext0=qtext0+gaprq+gapqq+qtext[i][s00:]
				rtext0=rtext0+gaprr+gapqr+rtext[i][newpposi:]
		

		return qstart0, qend0, qtext0, rtext0

	try:
		t=pd.read_csv(lastzname, skiprows=1 ,sep ='\t',header=None)
	except:
		return 0,0	
	
	
	frontfix=start-searchs

	t=t[(t[0]!='#name1')]

	l=len(t)
	
	if l<1:
		return 0,0
	if len(exons)<2:
		return 0,0
	
	qstart0=list(t[[1]].values.flatten())
	qend0=list(t[[2]].values.flatten())
	seq0=list(t[[10]].values.flatten())
	seq1=list(t[[11]].values.flatten())
	rstarto=list(t[[6]].values.flatten())
	rendo=list(t[[7]].values.flatten())
	strand0=list(t[[8]].values.flatten())
	chrom0=list(t[[5]].values.flatten())
	name0=list(t[[0]].values.flatten())
	score0=list(t[[12]].values.flatten())
	size0=list(t[[9]].values.flatten())
	nm=list(t[[13]].values.flatten())
	mm=list(t[[14]].values.flatten())
	size1=list(t[[15]].values.flatten())
	
	exonsfix=[searchs]+ exons[1:]
	exonefix=exone[:-1]+[searche]
	
	
	chrom1,rstart0,rend0,pieces=recoverinfor(chrom0,rstarto,rendo,size1)

	mrnasize, posi=findposi0(exons,exone)
	posi=[a+frontfix for a in posi] 
	
	index,gaps=findcombine(rstart0,rend0,qstart0,qend0, chrom0,strand0,5,seq0,seq1)
	#index=[[a] for a in xrange(len(qstart0))]


	
	extend=[0,0,[],[]]
	outputs=[]
	output=[]
	num=0
	leftnum=0
	if len(index)>0:

		for j,index1 in enumerate(index):
	      		
			if len(index1)==1:
				i=index1[0]
				if rend0[i]-rstart0[i]<=40:
					continue
					
				a0=[a for a in xrange(len(posi)) if posi[a]<=qstart0[i]+20]
				a1=[a for a in xrange(len(posi)) if posi[a]<qend0[i]-20]
				a2=[a for a in a1 if a not in a0]
				a2.sort()
				a2=[min(a2)-1]+a2		
				chrom10=str(chrom1[i])
	                
				absexonposi=cutrange([qstart0[i],qend0[i]],exonsfix, exonefix)
				absexonposi=[str(a) for a in  absexonposi]
				if qstart0[i] < 20 and searchs>0 :
					extend[0]=1
					
				if qend0[i] > size0[i]-20 and searche<lread:
					extend[1]=1
					
				if rstarto[i]<20 and pieces[i]>0:
					extend[2].append(str(chrom0[i]))
				
				if rendo[i]>size1[i]-20 and pieces[i]+size1[i]<piecesize[str(chrom1[i])]:
					extend[3].append(str(chrom0[i]))
						
					
				output.append([gene,reffile,chrom,int(start),int(end),len(exons),';'.join([str(a) for a in a2]), ';'.join([absexonposi[2*a] for a in xrange(len(absexonposi)/2)]),';'.join([absexonposi[2*a+1] for a in xrange(len(absexonposi)/2)]),qstart0[i]-frontfix,qend0[i]-frontfix,chrom10,rstart0[i],rend0[i],int(end)-int(start),mrnasize,rend0[i]-rstart0[i],float(qend0[i]-qstart0[i])*100.0/mrnasize,score0[i],nm[i],mm[i],float(nm[i])*100.0/(qend0[i]-qstart0[i]),strand,strand0[i],seq0[i],seq1[i]])
				leftnum= leftnum+1
					
				
			elif len(index1)>1:
				
				i=index1[0]			
				rstart10=[int(rstart0[a]) for a in index1]
				rend10=[int(rend0[a]) for a in index1]
				
				refposis=min(rstart10)
				refposie=max(rend10)
				a10=[]
				if refposie-refposis<=40:
					continue
					
				exonindex=[]
				s0=0
				for i in index1:
					
				    if qend0[i]-qstart0[i]<=40:
					continue
						
					a0=[a for a in xrange(len(posi)) if posi[a]<=qstart0[i]+20]
					a1=[a for a in xrange(len(posi)) if posi[a]<qend0[i]-20]
					a2=[a for a in a1 if a not in a0]
					a2.sort()			
					if a2 != []:
						a2=[min(a2)-1]+a2    					
						exonindex.extend(a2)
						s0=1
	                        	

				
				a2=list(set(exonindex))
				a2.sort()
	                    
                
				qstart10=[int(qstart0[a]) for a in index1]
				qend10=[int(qend0[a]) for a in index1]
				
				seq00=[seq0[a] for a in index1]
				seq01=[seq1[a] for a in index1]
	                
				qstart1,qend1,seq00,seq01=cordicombine(qstart10,qend10,rstart10,rend10,seq00,seq01,strand0[i])

				if ':' in seq00:
					with open('temp/{:s}/{:s}.fa'.format(reffile,gene) ,mode='r') as f:
						mseq0=''.join(f.read().split('\n')[1:])
					f.close()

					reposi=re.findall(r'\d+:\d+', seq01)
			
					for reposi0 in reposi:
						reposi0s=int(reposi0.split(':')[0])
						reposi0e=int(reposi0.split(':')[1])
						seq00=seq00.replace(reposi0,mseq0[reposi0s:reposi0e])

	                
				absexonposi=cutrange([x for y in [[qstart1[a],qend1[a]] for a in xrange(len(qstart1))] for x in y],exonsfix, exonefix)
					
	                	absexonposi=[str(a) for a in  absexonposi]
				i=index1[0]			
	               			
				chrom10=str(chrom1[index1[0]])
	                
				if min(qstart1) < 20  and searchs>0 :

					extend[0]=1
					
				if max(qend1) > size1[i]-20 and searche<lread:
					extend[1]=1
					
				if min([rstarto[a] for a in index1])<20 and pieces[i]>0:
					extend[2].append(chrom0[i])
				
				if max([rendo[a] for a in index1])>size1[i]-20 and pieces[i]+size1[i]<piecesize[str(chrom1[i])]: 
					
					extend[3].append(chrom0[i])
				if len(seq00.replace('-','')) != sum([int(absexonposi[2*k+1])-int(absexonposi[2*k]) for k in xrange(len(absexonposi)/2) ]):
					with open('errorlog', mode='a' ) as f:
						f.write(gene+'\n')
					f.close()

	
	                    
				output.append([gene,reffile,chrom,int(start),int(end),len(exons),';'.join([str(a) for a in a2]), ';'.join([absexonposi[2*a] for a in xrange(len(absexonposi)/2)]),';'.join([absexonposi[2*a+1] for a in xrange(len(absexonposi)/2)]),';'.join([str(a-frontfix) for a in qstart1]),';'.join([str(a-frontfix) for a in qend1]),chrom10,refposis,refposie,int(end)-int(start),mrnasize,refposie-refposis,float(refposie-refposis)*100.0/mrnasize,sum([score0[a] for a in index1]),sum([nm[a] for a in index1]),sum([mm[a] for a in index1]),float(sum([nm[a] for a in index1]))*100.0/(refposie-refposis),strand,strand0[i],seq00,seq01])
				    
	                    
				leftnum= leftnum+1
    
	if extend!=[0,0,[],[]] and cycle<6:
		return -1,extend
	       			
	allrecord=cl.defaultdict(list)

	for i,a in enumerate(output):

		allrecord['{:s}_{:s}_{:s}_{:s}_{:s}'.format(str(a[9]),str(a[10]), str(a[11]),str(a[12]),str(a[13]))].append(i)

	output0=[]
	for a in allrecord.values():
		output0.append(output[a[0]])
	
	leftnum=len(allrecord.values())
	output0=pd.DataFrame.from_records(output0)

	lock1.acquire()
	f = output0.to_csv("result/{:s}_retro.csv".format(reffile),sep=',', header=None, mode="a")
	outputs=pd.DataFrame.from_records(outputs)
	f = outputs.to_csv("result/{:s}_retrosp.csv".format(reffile),sep=',', header=None, mode="a")
	lock1.release()
	del t 
	del output
	del outputs
	
	gc.collect()

	return leftnum,num	



def runlastz(j,reffile=1,c0=0,chrom=0, gene=0,start=0,end=0,exons=0,exone=0,strand=0,lread=0,piecesize=0,extendloop=[],mode=1,cycle0=1):
		
	
	def extendfile(reffile,chrom,mRNAfile0,piecefile,extend,searchs,searche,exons,exone):
		
		
		#lock2.acquire()

			
		searchs0=searchs
		searche0=searche

		if extend[0]==1:
			searchs=max(0, searchs-1000)
		if extend[1]==1:
			searche=min(lread, searche+1000)
		
		if extend[0]+extend[1]>0:
			
			cm='python retrocm.py -f {:s} -c {:s}  -g {:s} -m 2 -s {:s} -e {:s} -x {:s} -y {:s}'.format(reffile.split('/')[-1],chrom, gene,'0'+ str(searchs), '0'+str(searche),'0'+ '#'.join([str(a) for a in exons]),'0'+ '#'.join([str(a) for a in exone]))
				
			os.system(cm)

		
		#lock2.release()
		
	
		return searchs, searche
	
	
	if type(j)==type([]):

		mode=1
	
		j,reffile,c0,chrom, gene,start,end,exons,exone,strand,lread,piecesize,extendloop,cycle0=j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8],j[9],j[10],j[11],j[12],j[0][1],j[13]			

	if mode==0:

		lock.acquire()
		
		if j not in runnedloop[c0]:
			runnedloop[c0]=runnedloop[c0]+[j]

		lock.release()
	
	mRNAfile0='temp/{:s}/{:s}.fa'.format(reffile,gene)
	
	piecefile='temp/{:s}/ref_{:s}.fa'.format(reffile,gene)
	
	outputfile='temp/{:s}/out_{:s}.csv'.format(reffile,gene)
	
	if extendloop !=[]:

		with open('runlog.txt',mode='a') as f:
			f.write('extend assembly range, gene:{:s}\n'.format(gene))
		f.close()	
		
		cycle=extendloop[1]
		searchs=extendloop[2]
		searche=extendloop[3]
			
	else:
		cycle=1
		searchs=max(0, start-500)
		searche=min(lread, end+500)
	
	while cycle>0 and cycle<7:
		
		if mode=='0':		
			with open('runlog.txt', mode='a') as f:
				f.write(chrom+' '+gene+' '+str(len(finishedloop[c0]))+'\n')
			f.close
		
		refsize=os.path.getsize(piecefile)/1900000000

		if refsize>0:
			cm='lastz_32 --format=general:name1,zstart1,end1,strand1,size1,name2,zstart2+,end2+,strand2,size1,text1,text2,score,nmatch,nmismatch,size2 --ambiguous=n --strand=both --identity=90 --ydrop=1000 --match=1,5 --exact=10  {:s} {:s} > {:s}'.format(mRNAfile0, piecefile, outputfile)

		else:
			cm='lastz --format=general:name1,zstart1,end1,strand1,size1,name2,zstart2+,end2+,strand2,size1,text1,text2,score,nmatch,nmismatch,size2 --ambiguous=n --strand=both --identity=90 --ydrop=1000 --match=1,5 --exact=10  {:s} {:s} > {:s}'.format(mRNAfile0, piecefile, outputfile)

		
		try:
			os.system(cm)
			s=0
		except:
			s=3
		
		while s>0:
			try:
				os.system(cm)
				s=0
			except:
				s=s-1
				with open('errorlog.txt', mode='a') as f:
					f.write('\n repeat lastz {:s} #{:d}\n'.format(gene,4-s))
				f.close
		
		
		num0,num1=readlastz(outputfile,reffile,c0,chrom,gene,start,end,exons,exone,strand,piecesize,searchs,searche,lread,max(cycle,cycle0))


		if num0>0:
	
			cycle=0
	        
			print 'Found retros in {:s} on chrom {:s}, gene {:s}, number of retro: #{:d}, assembled #{:d}\n'.format(reffile,chrom,gene,num0,num1)
	                
		elif num0==0:
		
			cycle=0
	                
			print 'Comfired no retron in {:s} on chrom {:s}, gene {:s}\n'.format(reffile,chrom,gene)
		
		else:

			print 'Extend retro #{:d}  in {:s} on chrom {:s}, gene {:s}\n'.format(cycle,reffile,chrom,gene)		
			if cycle>6:
				break	
			cycle=cycle+1
			
			searchs,searche=extendfile(reffile,chrom, mRNAfile0,piecefile,num1,searchs, searche,exons,exone)
			
			if num1[2] !=[] or  num1[3] !=[]:
				
				try:
					os.remove(outputfile)
				except:
					print ('unable remove {:s}\n'.format(outputfile)) 
				
				extendloop=[gene,cycle,searchs, searche]
				
				if mode==0:	
					lock.acquire()  
					if j not in finishedloop[c0]:
						finishedloop[c0]=finishedloop[c0]+[j]
						print len(finishedloop[c0])
					lock.release()
				
				return [[num1,extendloop],j,reffile,c0,chrom, gene,start,end,exons,exone,strand,lread,piecesize]
				
			
	with open('runlog.txt',mode='a') as f:
		f.write('\n\tfinised {:s} {:s} #{:d}'.format(reffile,chrom,j)+'\n')
	f.close	
			

	try:
		os.remove(mRNAfile0)
	except:
		print ('unable remove {:s}\n'.format(mRNAfile0)) 
		
	
	try:
		os.remove(outputfile)
	except:
		print ('unable remove {:s}\n'.format(outputfile)) 
			

	try:
		os.remove(piecefile)
	except:
		pass


	if mode==0:	
		lock.acquire()  
		if j not in finishedloop[c0]:
			finishedloop[c0]=finishedloop[c0]+[j]
			print len(finishedloop[c0])
		lock.release()

		
	return 0
		
def runhit(reffiles,allchromfiles):	
	
	def checkunfinished(reffiles0):
		
		hitstart={}
		
		try:
			allfiles=os.listdir(savefolder+'result/'+reffiles0+'/')
		except:
			allfiles=[]
		finishedchroms=['_'.join(allfiles[a].split('_')[:-1]) for a in xrange(len(allfiles)) if os.path.isfile(os.path.join(savefolder+'result/'+reffiles0+'/',allfiles[a])) and re.search('exonout.csv$',allfiles[a])!=None]
		
		allfinishedchroms=[]
		exonnumbers={'chr20.fa': 12801, 'chr3.fa': 32480, 'chr4.fa': 22157, 'chrY.fa': 1409, 'chr6.fa': 25655, 'chr10.fa': 21308, 'chr19.fa': 30878, 'chr12.fa': 31142, 'chr1.fa': 52758, 'chr21.fa': 6434, 'chr8.fa': 19754, 'chr18.fa': 9441, 'chr17.fa': 31497, 'chr15.fa': 19934, 'chr7.fa': 25603, 'chrX.fa': 16683, 'chrM.fa': 0, 'chr13.fa': 9696, 'chr11.fa': 30031, 'chr16.fa': 23877, 'chr14.fa': 17782, 'chr9.fa': 20347, 'chr2.fa': 41969, 'chr5.fa': 24750, 'chr22.fa': 11412}	
	
		for finishedchrom0 in finishedchroms:


			try:
				t=pd.read_csv('result/{:s}/{:s}_exonout.csv'.format(reffiles0,finishedchrom0),header=None,sep='\t')
				
				last=list(t[1])[-1]
				del t
			except:
				last=[]
			
			if last!=[]:
				
				genefile=exonfolder+'{:s}_gene.csv'.format(finishedchrom0.split('.fa')[0])
			
				exonlist=pd.read_csv(genefile,skiprows=0,sep='\t',header=None)

				t1=list(exonlist[2])
				last0=(t1.index(last)+1)/10000

				if last0==len(t1)/10000:
					
					allfinishedchroms.append(finishedchrom0)			

				else:

					hitstart[finishedchrom0]=last0*10000
					
					
		
		
		return allfinishedchroms+['chrM.fa'], hitstart
	
	
	record=[]
	for reffiles0 in reffiles:
		
		if unfinished==1:
			finishedchrom,hitstart=checkunfinished(reffiles0)
			chromfiles=[a for a in allchromfiles if a not in finishedchrom]
		
		else:
			chromfiles=allchromfiles
			hitstart={}
			
		if chromfiles==[]:
			continue
		
		print 'running assembly {:s} with choromosome {:s}\n'.format(reffiles0,str(chromfiles))
		
		with open('runlog.txt',mode='a') as f:
			f.write('\n\tstart '+reffiles0+'\n')
		f.close
		
		m0=1
		p=mul.Pool(processes=m0)
		record0=[]

		for i,chrom in enumerate(chromfiles):
			
			hitstart0=0
			if chrom in hitstart.keys():
				
				hitstart0=hitstart[chrom]
			
			time.sleep(i-i/m0*m0)
			print '\nrunning reference file: {:s}\n'.format(chrom)
			record0.append([])
			record0[i]=p.apply_async(hit,(i,reffiles0,chrom,hitstart0))
		p.close()
		p.join()
		p=mul.Pool(processes=m0)
		s=[3 for a in xrange(len(chromfiles))]
		while min(s)>0:
			
			for i,chrom in enumerate(chromfiles):
				
				if s[i]==0:
					continue	
				try:
					record0[i]=record0[i].get()			
					s[i]=0
				except:
					record0[i]=p.apply_async(hit,(i,reffiles0,chrom))
					s[i]=s[i]-1
					with open('errorlog.txt', mode='a') as f:
						f.write('\n repeat hit{:s} {:s} #{:d}\n'.format(reffiles0,chrom,4-s[i]))
					f.close
		p.close()
		p.join()
		record.append(record0)

	return record
	
		
	
	
	


def gu(reffiles,chrom):
	time.sleep(2)
        try:
		a=runned.index(os.getpid())
        except:
		a=-1

        return a


def eachchromtranstogene(reffile,chrom,trans,genes,exonnotes,transname,exons,exone,strand):

	file0='result/{:s}/{:s}_exonout.csv'.format(reffile,chrom)
	
	if  '{:s}_exonout.csv'.format(chrom) not in os.listdir('result/{:s}/'.format(reffile)):
		return 0 
	
	list0=pd.read_csv(file0,sep='\t',header=None)
	trans0=list(list0[1])
	posi=list(list0[4])

	find=fasterfind([x.split('.')[0] for x in trans0],trans)
	genelist=[genes[i] for i in find]
	
	indextrans=[transname.index(x) for x in trans0]
	chromstrand={genelist[i]:strand[indextrans[i]] for i in xrange(len(genelist))}
	transposis={genelist[i]:';'.join(sorted(exons[indextrans[i]].split(',')[:-1]+exone[indextrans[i]].split(',')[:-1])) for i in xrange(len(genelist))}


	common=cl.defaultdict(list)
	closeexon=cl.defaultdict(list)	
	for i,x in enumerate(genelist):
	
	 	if type(posi[i])!=type('nan'):
			continue
		
		posi0=posi[i].split(';')
		pieces=[a.split(':')[0] for a in posi0]
		posis=[a.split(':')[1] for a in posi0]
		
		if len(posis)>1000:
			
			for j,piece in enumerate(pieces):
				common[x+'_'+piece].append(int(posis[j]))
			continue

		for j,piece in enumerate(pieces):
			closeexon[x+'_'+piece].append(int(posis[j]))
		
		
	
	allhits=cl.defaultdict(list)
	for piece in closeexon.keys():

		
		hitpiece=piece.split('_')[1]
		posis=sorted(closeexon[piece])
		
		if piece in common.keys():
			
			posis2=sorted(common[piece])
			
			for posi0 in posis2:
				
				if min([abs(a-posi0) for a in posis])<30:
					
					posis.append(posi0)	

		posis.sort()			
		
		distances=[posis[i+1]-posis[i] for i in xrange(len(posis)-1)]
		

		close=sorted(list(set([posis[i]+distances[i]/2 for i in xrange(len(distances)) if distances[i]<30])))
		
		if close==[]:
			continue

		allhits[piece.split('_')[0]].extend([[hitpiece,close0] for close0 in close])
	
	for hitposis in allhits.keys():
		
		
		allhits[hitposis]=';'.join([str(a[0])+':'+str(a[1]) for a in sorted(allhits[hitposis],key=lambda x:x[0])])
		
	out=[[genename, allhits[genename],chromstrand[genename], exonnotes[genename]] if genename in exonnotes.keys() else [genename,allhits[genename],chromstrand[genename], transposis[genename]]  for genename in allhits.keys()]
	
	write=pd.DataFrame.from_records(out).to_csv('result/{:s}/{:s}_geneout.csv'.format(reffile,chrom),sep=',',header=None,mode="w")
		
	
	return 0



def transtogene(reffiles,chromfiles):
		
		transfile='transinfor.csv'
		genesfile='geneinfor.csv'
		
		
		data=pd.read_csv(transfile,sep='\t',header=None)
		trans=list(data[1])
		genes=list(data[2])		

		
		data1=pd.read_csv(genesfile,sep='\t',header=None)
		genename=list(data1[1])
		exonposis=list(data1[2])
		
		exonnotes={}
		for i,genename0 in enumerate(genename):
			exonnotes[genename0]=exonposis[i]
			
		
		del data,data1,genename, exonposis
		
		list1=pd.read_csv('knownGenefix.txt',sep='\t',header=None)
		transname=list(list1[12].values.flatten())
			
		exons=list(list1[9].values.flatten())
		exone=list(list1[10].values.flatten())
		
		strand=list(list1[3].values.flatten())
		del list1
		
		
		for reffile in reffiles:
			
			m0=15
			p=mul.Pool(processes=m0)
			transtogene_record=[]	
			for i,chrom in enumerate(chromfiles):	
				oldname='result/{:s}/{:s}_geneout.csv'.format(reffile,chrom)
				newname='result/{:s}/{:s}_tranout.csv'.format(reffile,chrom)
				try:
					os.system('mv {:s} {:s}'.format(oldname,newname))
				except:
					pass
				transtogene_record.append(p.apply_async(eachchromtranstogene,(reffile,chrom,trans,genes,exonnotes,transname,exons,exone,strand)))
				

			p.close()
			p.join()

			for record0 in transtogene_record:
				print record0.get()
		

		









		
def psize(reffile):
	
	with open(wkfolder+reffile,mode='r') as f:
		reado=f.read()
	f.close()


	read=reado.split('>')[1:]

	del reado
	gc.collect()

	psize={}
	for read0 in read:
		cut=read0.split('\n')
		title=cut[0].split(' ')[0]
		size=sum([len(x) for x in cut[1:]])
		psize[str(title)]=size+1-1
		del read0

	del read		

	gc.collect()

	return psize
	

def extendref(reffile,refseq,reftitle,gene,extend):
	if extend==[]:
		return 0
	if refseq==[]:
		with open(wkfolder+reffile,mode='r') as f:
			read0=f.read().split('>')[1:]
		f.close()
		reftitle=[a.split('\n')[0].split(' ')[0] for a in read0]
		refseq=[''.join(a.split('\n')[1:]) for a in read0]
		del read0



		
	piecefile='temp/{:s}/ref_{:s}.fa'.format(reffile,gene)	

	
	with open(piecefile, mode='r') as f:
		read=f.read().split('>')[1:]
	f.close()	
			
	title=[a.split('\n')[0] for a in read]
		
	a2=[title.index(a) for a in extend[2]]
	a3=[title.index(a) for a in extend[3]]
		
	write=''
	for x in xrange(len(title)):
			
		if x not in a2 and x not in a3:
			write=write+'>'+read[x]+'\n'
			continue
			
		a=title[x]
		pindex=a.split('_')[0]
		pstart=int(re.search(r'(?<=_)\d+(?=(:|;))',a).group())
		pend=int(re.search(r'(?<=;)\d+',a).group())
			
		refseq0=refseq[reftitle.index(pindex)]+''
			
		if x in a2:
				
			pstart=max(pstart-5000, 0)
			
		if x in a3:
				
			pend=min(pend+5000, len(refseq0))
				
		write=write+'>{:s}_{:d}_{:d};{:d}\n'.format(pindex,x,pstart,pend)+refseq0[pstart:pend]+'\n'
			
		del refseq0
		
		
	with open(piecefile,mode='w') as f:
		f.write(write)
	f.close()
	
	del write
	

	
def makeup(reffile, left):

	global runnedloop, finishedloop	
	
	runnedloop=[0 for x in xrange(100)]
	finishedloop=[0 for x in xrange(100)]

	time0=0
	makeup=[x for x in range(len(left))]

		
	while len(makeup)>0 and time0<3:

		time0=time0+1
		makeup0=[]
		allcm=''
		if type(left)!=type([]):
			left=left.get()
		
		for i in xrange(len(left)):
		
			if type(left[i])!=type([]) and  type(left[i])!=type(0):			
				try:
					left[i]=left[i].get()
				except:
					print left[i]
			if type(left[i])==type([]):
				
				if left[i][0]!=[]:
				
					cms='0'+'#'.join(left[i][0][0][2])
					cme='0'+'#'.join(left[i][0][0][3])
					
					cm='python retrocm.py -w {:s}  -f {:s} -g {:s} -m 0 -s {:s} -e {:s}'.format(wkfolder,reffile, left[i][5],cms,cme)
					allcm=allcm+left[i][5]+'!'+cms+'!'+cme+'\n'
					#extendref(left0[2],refseq,reftitle,left0[5],left0[0][0])

				makeup0.append(i)

		with open(reffile+'_transcm',mode='w') as f:
			f.write(allcm)
		f.close()

		cm='python retrocm.py -w {:s}  -f {:s}  -m 0'.format(wkfolder,reffile)
		os.system(cm)
		gc.collect()

		m0=15
		p=mul.Pool(processes=m0)
		left=p.map_async(runlastz,[left[i]+[time0+3] for i in makeup0])
			
		p.close()
		p.join()
			

	try:
		checkmem('makeup')
	except:
		pass

	cm='python retrocm.py -w {:s} -f {:s} -m 1'.format(wkfolder,reffile)
	os.system(cm)
	
	checkmem('fix')
	
	gc.collect()



		
	




def runreffile(reffile,chromfiles,allpsize):
				
	

	def readhitfile(genefile,reffile,chrom):
		#lock2.acquire()

		try:
			geneinfo=pd.read_csv(genefile,sep=',',skiprows=0,header=None)
		except:
			return  [],[],[],[],[],[],[],[]
		
		print 'start analyzing reference # {:s} with chrom {:s}, has candidates {:d}\n'.format(reffile,chrom,len(geneinfo))
		
		genesort=list(geneinfo[0].values.flatten())
		genename=list(geneinfo[1].values.flatten())
		hitposi=list(geneinfo[2].values.flatten())
		strand=list(geneinfo[3].values.flatten())
		exonse=[a.split(';') for a in list(geneinfo[4].values.flatten())]
		

		exons=[[int(b[2*k]) for k in xrange(len(b)/2)] for b in exonse]
		exone=[[int(b[2*k+1]) for k in xrange(len(b)/2)] for b in exonse]
		
		
		start=[min(a) for a in exons]
		end=[max(a) for a in exone]


		hitpiece=[[b for b in re.findall(r'(\d+|\w+|_)+(?=:)',a)] for a in hitposi]
		hitposi=[[int(b) for b in re.findall(r'(?<=:)[0-9]+',a)]  for a in hitposi]

		


		return genename,hitpiece, hitposi,start,end,exons,exone,strand
		
	
	def eachchrom(c0,reffile,chrom,p,left,finishedindex,piecesize):
		global runnedloop, finishedloop
		
		with open('runlog.txt',mode='a') as f:
			f.write('\nanalyze '+reffile+' '+chrom+'\n')
	 	f.close
		
		genefile='result/{:s}/{:s}_geneout.csv'.format(reffile,chrom)
		
		exonfile='result/{:s}/{:s}_exonout.csv'.format(reffile,chrom)
	
		cm='python retrocm.py -w {:s} -f {:s} -c {:s} -m 3'.format(wkfolder, reffile, chrom)
		gc.collect()

		print cm
		os.system(cm)

		#checkmem('os')
		alllread={'chr20.fa': 64444167, 'chr3.fa': 198295559, 'chr4.fa': 190214555, 'chrY.fa': 57227415, 'chr6.fa': 170805979, 'chr10.fa': 133797422, 'chr19.fa': 58617616, 'chr12.fa': 133275309, 'chr1.fa': 248956422, 'chr21.fa': 46709983, 'chr8.fa': 145138636, 'chr18.fa': 80373285, 'chr17.fa': 83257441, 'chr15.fa': 101991189, 'chr7.fa': 159345973, 'chrX.fa': 156040895, 'chrM.fa': 16569, 'chr13.fa': 114364328, 'chr11.fa': 135086622, 'chr16.fa': 90338345, 'chr14.fa': 107043718, 'chr9.fa': 138394717, 'chr2.fa': 242193529, 'chr5.fa': 181538259, 'chr22.fa': 50818468}
		lread=alllread[chrom]
		genename,hitpiece, hitposi,start,end,exons,exone,strand=readhitfile(genefile,reffile,chrom)

		gc.collect()

		try:
			checkmem('read')
		except:
			pass


		record12=[]
		#m0=16
		#p=mul.Pool(processes=m0) 
		for j, gene in enumerate(genename):
			#j=genename.index('ENST00000399796.6')
			#gene='ENST00000399796.6'
			#if j >=50:
			#continue
			
			#runlastz(j,reffile,c0,chrom, gene,start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize)
			#runlastz(j,reffile,c0,chrom, gene,start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize)
			record12.append(p.apply_async(runlastz,(j,reffile,c0,chrom, gene,start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize,[],0)))

			
		#p.close()
		#p.join()
		#for r in record12:
			#r.get()

		time.sleep(10)
		success=[]
		while len(runnedloop[c0])<len(genename)-2:
			with open ('current.out',mode='w') as f:

				f.write(str(c0)+' '+chrom+'\n'+'running:'+str([a for a in runnedloop[c0] if a not in finishedloop[c0]])+ '\nfinished:'+str(finishedloop[c0]))
			f.close()

			time.sleep(30)

		finishedindex.append([a for a in finishedloop[c0]])
		left.extend([record12[j] for j in xrange(len(record12)) if j not in finishedindex[-1]])

		getted=[a for a in finishedindex[-1]] 
		gettednew=[]
		for j in getted:
 				
			record12[j]=record12[j].get()
			
			if type(record12[j])==type([0]):
				print 'extend', j
				gettednew.append(j) 

		left.extend([record12[j] for j in gettednew])
		
		getted=[]
		#getted=[a for a in finishedindex[-1]]			
		if  len(getted) >0:

			with open(wkfolder+reffile,mode='r') as f:
				read0=f.read().split('>')[1:]
			f.close()
		
			reftitle0=[a.split('\n')[0].split(' ')[0] for a in read0]
			refseq0=[''.join(a.split('\n')[1:]) for a in read0]
			del read0
			gc.collect()

		while len(getted) >0:

			gettednew=[]
			for j in getted:
	
				record12[j]=record12[j].get()
					
				
				if type(record12[j])==type([0]):
					print 'extend', j
					gettednew.append(j)
					if record12[j]!=[]:
						cm='python retrocm.py -w {:s} -f {:s} -g {:s} -m 0 -s {:s} -e {:s}'.format(wkfolder,reffile, genename[j], ','.join(record12[j][0][2]),','.join(record12[j][0][3]))
						os.system(cm)

						#extendref(reffile,refseq0,reftitle0,genename[j],record12[j][0])

			for j in gettednew:
				
				if record12[j]!=[]:
					
					record12[j]=p.apply_async(runlastz,(j,reffile,c0,chrom, genename[j],start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize,record12[j][1]))
				else:
		
					record12[j]=p.apply_async(runlastz,(j,reffile,c0,chrom, genename[j],start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize))
					
					
			getted=gettednew	
				


				#break 
				#try:
					#r00=record12[j].get()
					#while type(r00)==type([0]):
						#extendref(reffile,genename[j],r00)
						#record12[j]=p.apply_async(runlastz,(j,reffile,c0,chrom, genename[j],start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize))
						#r00=record12[j].get()
					#success=1
				#except:
									
					#with open('errorlog.txt',mode='a') as f:
						#f.write(str(repeat) +' '+gene+'\n')
					#f.close
					#repeat=repeat-1
					#record12[j]=p.apply_async(runlastz,(j,reffile,c0,chrom, genename[j],start[j],end[j],exons[j],exone[j],strand[j],lread,piecesize))
				
		map(deepdel,[record12,genename,hitpiece, hitposi,start,end,exons,exone,strand,lread,piecesize])

		return finishedindex,left

		
	retrofile='./result/{:s}_retro.csv'.format(reffile)
	retrospfile='./result/{:s}_retrosp.csv'.format(reffile)

	try:
		os.system('mv {:s} {:s}_bak'.format(retrofile,retrofile))
	except:
		pass

	try:
		os.system('mv {:s} {:s}_bak'.format(retrospfile,retrospfile))
	except:
		pass


	header= ','.join(['sortnumber','genename','assembly','chr','gene_start_on_chr','gene_end_on_chr','involve_exons','involve_exons_start_posi','involve_exons_end_posi','retro_start_on_mrna','retro_end_on_mrna', 'retro_on_assemble_piece_title', 'retro_on_assemble_piece_start', 'retro_on_assemble_piece_end','gene_size','alignment_size', 'alignment_score', 'retro_strand','alignment_strand','text1','text2'])


	f=open(retrofile,mode='w')
	f.write(header+'\n')
	f.close()
	
	f=open(retrospfile,mode='w')
	f.write(header+'\n')
	f.close()

	if unfinished==1:
		try:
			allfiles=os.listdir(savefolder+'result/'+reffile+'/')
			finishfiles=[allfiles[a].split('_')[0] for a in xrange(len(allfiles)) if os.path.isfile(os.path.join(savefolder+'result/'+reffile+'/',allfiles[a])) and re.search('geneout.csv$',allfiles[a])!=None]
			finishfiles2=[allfiles[a].split('_')[0] for a in xrange(len(allfiles)) if os.path.isfile(os.path.join(savefolder+'result/'+reffile+'/',allfiles[a])) and re.search('exonout.csv$',allfiles[a])!=None]
			finishfiles=[a for a in finishfiles if a in finishfiles2]
		except:
			finishfiles=[]
		chromfiles=finishfiles	

	
	global finishedloop,runnedloop
	finishedloop=manager.list([[] for a in chromfiles])
	runnedloop=manager.list([[] for a in chromfiles])
	

	left=[]	
	m0=15
	p=mul.Pool(processes=m0) 			
	finishedindex=[]
	for i,chrom in enumerate(chromfiles):	
		
		finishedindex,left=eachchrom(i,reffile,chrom,p,left,finishedindex,allpsize)

	p.close()
	p.join()

	map(deepdel,[finishedloop,runnedloop])	
	makeup(reffile,left)	

	deepdel(finishedindex)

	gc.collect()	
	try:
		checkmem('reffile')
	except:
		pass
		

def runanalyze(reffiles,chromfiles):
	
	for reffile in reffiles:

		try:
			checkmem('start')
		except:
			pass

		allpsize=psize(reffile)

		runreffile(reffile,chromfiles,allpsize)

		gc.collect()

		try:
			checkmem('end')
		except:
			pass
			
	return 0

	
	




















def needleallfiles(i,tempfolder,alltrans,chroms,refseq,try0=0):
	
	
	print alltrans[0]	
	mstart0=[[tran[2*a] for a in xrange(len(tran)/2)] for tran in alltrans]
	mend0=[[tran[2*a+1] for a in xrange(len(tran)/2)] for tran in alltrans]		
	
	mfile=tempfolder+'mrnax_{:d}.fa'.format(i)
	
	rfile=tempfolder+'refx_{:d}.fa'.format(i)
	
	ofile=tempfolder+'multi_{:d}.out'.format(i)
	
	mrnasize=[]
	errorindex=[]

	with open(mfile, mode='w') as f:

		for j in xrange(0,len(mend0)):

			if sum(mend0[j])-sum(mstart0[j])==0:
				errorindex.append(j)
				mrnasize.append(0)
				continue
			
				
			mrnasize.append(sum(mend0[j])-sum(mstart0[j]))
			
			
			f.write('>m{:d}\n'.format(j)+''.join([chroms[mstart0[j][a]:mend0[j][a]] for a in xrange(len(mend0[j])) ])+'\n')
				
	f.close()
	

	with open(rfile, mode='w') as f:
		f.write('>s0\n'+refseq.replace('-',''))
	f.close()
	
	return mrnasize


def writecomparefiles(i,tempfolder,alltrans,chroms,refseq,highscore):
	
	mstart0=[[tran[2*a] for a in xrange(len(tran)/2)] for tran in alltrans]
	
	mend0=[[tran[2*a+1] for a in xrange(len(tran)/2)] for tran in alltrans]
	
	
	with open(tempfolder+'mrna_{:d}.fa'.format(i), mode='w') as f:
		
		f.write('>s1\n'+''.join([chroms[mstart0[highscore][a]:mend0[highscore][a]] for a in xrange(len(mend0[highscore]))]))		
	f.close()
	
	with open(tempfolder+'fullgene_{:d}.fa'.format(i), mode='w') as f:
		
		if mstart0[highscore]!=[]:
			
			f.write('>s1\n'+chroms[max(0,min(mstart0[highscore])):min(max(mend0[highscore]),len(chroms))])
			
		else:
			
			f.write('>s1\n'+'')
				
	f.close()
	
	
	try:
		os.remove(tempfolder+'mrnax_{:d}.fa'.format(i))
	except:
		pass
	

		
	
	
	
	
	


def findbest(i,tempfolder,alltrans,chrom0,refseq,mrnasize,try0=0):
		
		
	def readmulti(file0,l):
		print file0		
		with open(file0) as f:
			
			read=f.read().split('\n')[0:l]				
		f.close
		read=[a for a in read if a[0]=='s']		
			
		score=[float(re.search(r'\(\d+.*\d*\)', a).group().strip()[1:-1]) for a in read]
		
		return score

	mstart0=[[tran[2*a] for a in xrange(len(tran)/2)] for tran in alltrans]

	mend0=[[tran[2*a+1] for a in xrange(len(tran)/2)] for tran in alltrans]	
	
	mfile=tempfolder+'mrnax_{:d}.fa'.format(i)
	
	rfile=tempfolder+'refx_{:d}.fa'.format(i)
	
	ofile=tempfolder+'multi_{:d}.out'.format(i)
	
	errorindex=[a for a in xrange(len(mrnasize)) if mrnasize[a]==0]
	
	cm='needleall {:s} {:s} -snucleotide2  -gapopen 6 -gapextend 1  {:s}'.format(rfile,mfile,ofile)
	
	os.system(cm)
	
	if try0!=0:
		try:
			scores=readmulti(ofile,len(mend0))
			tempindex=0
			newscore=[]
			for a in xrange(len(mend0)):
				if a in errorindex:
					newscore.append(0.0)
				else:
					newscore.append(scores[tempindex])
					tempindex=tempindex+1
			scores=newscore	

			
		except:
			return -1,[0 for s in mrnasize],[]
	else:
		
		scores=readmulti(ofile,len(mend0)-len(errorindex))
		newscore=[]
		tempindex=0
		for a in xrange(len(mend0)):
			if a in errorindex:
				newscore.append(0.0)
			else:
				newscore.append(scores[tempindex])
				tempindex=tempindex+1
		scores=newscore

	mhighscore=max(scores[1:])
	 	
		
	if mhighscore>scores[0]-50:
		highscores=sorted([a for a in xrange(0,len(scores)-1) if scores[a+1]==mhighscore], key=lambda x: mrnasize[x+1])
		highscore=highscores[0]+1
	else:
		highscores=sorted([a for a in xrange(0,len(scores)-1) if scores[a+1]==mhighscore], key=lambda x: mrnasize[x+1])		
		highscore=0
	print i,sorted(scores[1:]), scores[0],highscore-1	
	
	try:
		os.remove(tempfolder+'multi_{:d}.out'.format(i))
	except:
		pass
	

	return highscore-1, scores,highscores




def combinemrna(file0,outputfile):


	def findoverlapi(start,end,index=''):
		
		l=len(start)
		allposi=start+end
		
		if index=='':
			index=range(l)
		key0=sorted(range(2*l), key=lambda i: allposi[i])
		c=0
		s=0
		start0,end0,index0,index1=[],[],[],[]
		for x in key0:
			
			if x < l:
				c=c+1
				index1.append(index[x])
				if c==1:
					s=start[x]
			else:
				c=c-1
				
			if c==0:
				start0.append(s)
				end0.append(end[x-l])
				index0.append(index1)
				index1=[]
				
		return start0, end0, index0
		
	
	def findpara(index,rstart,rend,score):
		
		
		rstart0,rend0,index0=findoverlapi(rstart,rend,index)
		
		index1,domi=[],[]
		for index10 in index0:
			
			index10=sorted(index10, key=lambda x: score[x],reverse=True)
			domi.append(index10[0])
			
			index1.append(index10)		
		
		return index1,domi
	

		
	def readinfor(file0,names='all'):
		
		t=pd.read_csv(file0, skiprows=0 ,sep ='\t',header=None)
		
		name=list(t[[12]].values.flatten())
		exons=list(t[[9]].values.flatten())
		exone=list(t[[10]].values.flatten())
		
		exons0=[]
		exone0=[]
		
		if names=='all':
			names=name
			
		for name0 in names:
			if name0 in name:
				i=name.index(name0)
				exons0.append([int(a) for a in exons[i].split(',')[:-1]])
				exone0.append([int(a) for a in exone[i].split(',')[:-1]])

		return exons0,exone0
	

	
	
	def combineretro(exonstart,exonend,name0,ref0,chrom0,strand1,piecei0,rstart0,rend0,gstart0,gend0,score0,seq0,seq1):

		def fixnega(qrange):
			
			exons=[qrange[2*k+1]-qrange[2*k] for k in xrange(len(qrange)/2)]
			trons=[qrange[k+1]-qrange[k] for k in xrange(len(qrange)-1)]
			
			shifts=0
			exonsize0=0
			exonsize=[]
			mark0=[]
			mark=[]
			
			for ioe,tron in enumerate(trons):
				
				if ioe-ioe/2*2!=0:
					
					if tron<0:
						
						shifte=ioe+1
						shifts=ioe+1
						
						for y in xrange(ioe/2):
							if qrange[ioe-y*2]<qrange[shifte]:
								shifts=ioe-y*2+1			
						
						mark1=range(shifts,shifte)
						mark0=[x for x in mark0 if x not in mark1]
						
						mark.append(mark0)
						mark.append(mark1)
						
						exonsize1=sum([exons[x/2] for x in mark1])/2
						exonsize.append(exonsize0-exonsize1)
						exonsize.append(exonsize1)			
						
						mark0=[]
						exonsize0=0	
										
						
						if tron<0:
							shifts=qrange[ioe]
						else:
							shifts=-1
				else:
					exonsize0=exonsize0+tron
					mark0.append(ioe)
					mark0.append(ioe+1)
					
				exonsize.append(exonsize0)
				mark.append(mark0)
				mark0=[]
				exonsize0=0
						
			clean0=[a for a in xrange(len(exonsize)) if exonsize[a]>20]
			mark=[x for a in clean0 for x in mark[a]]
			
			qrange0=[qrange[a] for a in mark]
			
			markposi=[[sum(exons[:x]),sum(exons[:x+1])] for x in xrange(len(exons)) if x not in list(set([a/2 for a in mark]))]
			
			return qrange0,markposi


				
		
		
		
		group=cl.defaultdict(list)
			
		for i in xrange(len(exonstart)):
			
			group[str(piecei0[i])+'&'+str(ref0[i])+'&'+chrom0[i]+'&'+strand1[i]+'&'+name0[i]].append(i)
			
			
			
		mstart=[]
		mend=[]
		indexname=[]
		indexsort=[]
		megarefs=[]
		megarefe=[]
		megagstart0=[]
		megagend0=[]
		megam=[]
		refchr=[]
		refseq=[]
		mrnaseq=[]
		addnames=[]
		
		sort00=0
		
		for k in group.keys():
			
			
			i=group[k]
			
			megarefs0,megarefe0,indexeachpiece=findoverlapi([rstart0[a] for a in i],[rend0[a] for a in i],i)
			
			megarefs.extend(megarefs0)
			megarefe.extend(megarefe0)
			indexsort.extend(indexeachpiece)
			
			refchr.extend([k.split('&')[0] for a in xrange(len(indexeachpiece))])
				
			strandx=k.split('&')[-2]
			
			for j,indexp0 in enumerate(indexeachpiece):
				indexp0=sorted(indexp0,key=lambda x:score0[x])
				indexname.append(name0[indexp0[0]])
				
				
			
				rrange=[[rstart0[indexp0[0]],rend0[indexp0[0]]]]

				qrange=[[a  for i in xrange(len(exonstart[indexp0[0]])) for a in [exonstart[indexp0[0]][i],exonend[indexp0[0]][i]]]]
				
				qpiece=[seq0[indexp0[0]].replace('-','')]
				rpiece=[seq1[indexp0[0]].replace('-','')]
				range0=[rstart0[indexp0[0]],rend0[indexp0[0]]]
				addnames0=[name0[indexp0[0]]]
				for a in indexp0[1:]:
					print k,addnames0,'retrorange=',rrange,'exonrange=',qrange,'\n\n','nextpiece','retrorange' ,[rstart0[a],rend0[a]],'generange', exonstart[a],exonend[a]

					addrrange=cutabsrange(rstart0[a],rend0[a],[0]+range0+[rend0[a]+1])

					print 'add',name0[a],'retrorange=',addrrange	
					if addrrange!=[]:
						addnames0.append(name0[a])
					else:
						continue
					if strandx=='+':
						addrrange0=[findcodi(seq1[a],[addrrange[k]-k+k/2*2-rstart0[a]])[0]+k-k/2*2 for k in xrange(len(addrrange))]						
					else:
						addrrange0=sorted([findcodi(seq1[a],[rend0[a]-addrrange[k]-k+k/2*2])[0]+k-k/2*2 for k in xrange(len(addrrange))])  

					addqrange0=[k-seq0[a][:k].count('-') for k in addrrange0]
				
					addqrange0=[[addqrange0[2*k],addqrange0[2*k+1]] for k in xrange(len(addqrange0)/2)]

					addqrange=[cutrange(k,exonstart[a],exonend[a]) for k in addqrange0]
					
					addrrange0=[[addrrange0[2*k],addrrange0[2*k+1]] for k in xrange(len(addrrange0)/2)]
					print 'exonrange=',addqrange

					if addqrange[-1][0] < qrange[-1][-1] and sum([addqrange[-1][2*x+1]-addqrange[-1][2*x] for x in xrange(len(addqrange[-1])/2)])<100:
						continue

					qrange.extend(addqrange)
					
					qpiece.extend([seq0[a][k[0]:k[1]].replace('-','') for k in addrrange0])					
					
					rpiece.extend([seq1[a][k[0]:k[1]].replace('-','') for k in addrrange0])
					
					range0=cleansame(sorted(addrrange+range0))
					print 'added', rrange,qrange,'\n'
					addrrange=[[addrrange[2*k],addrrange[2*k+1]] for k in xrange(len(addrrange)/2)]
					
					rrange.extend(addrrange)
					
				if strandx=='+':	
					key0=sorted(range(len(rrange)), key=lambda x: rrange[x][0])
				else:
					key0=sorted(range(len(rrange)), key=lambda x: rrange[x][0],reverse=True)

				rrange=range0
				
				piecesize=[sum([y[2*x+1]-y[2*x] for x in xrange(len(y)/2)]) for y in qrange]
				error=[x for x in xrange(len(key0)-1) if qrange[key0[x+1]][0]-qrange[key0[x]][-1]<0]
				errorindex=[x for x in xrange(len(piecesize)) if piecesize[x]<100 and (x-1 in error or x in error)]

					
				qrange=cleansame([a for x in key0 if x not in errorindex for a in qrange[x]])

				for x in errorindex:
					qpiece[x]=''.join(['-' for m in xrange(piecesize[x])])
				
				qrange,mark=fixnega(qrange)

				print '\n','endcombine','finalretrorange',rrange,'finalgenerange',qrange,'\n','finalexonsize/intronsize',[qrange[k+1]-qrange[k] for k in xrange(len(qrange)-1)]
											
				
				mrnaseq00=''.join([qpiece[x] for x in key0])

				for x in mark:
					mrnaseq00=mrnaseq00[:x[0]]+''.join(['-' for m in xrange(x[1]-x[0])])+mrnaseq00[x[1]:]

				
				refseq.append(''.join([rpiece[x] for x in key0]))
				mrnaseq.append(mrnaseq00)
				
				megam.append(qrange)
				addnames.append(addnames0)				
				
		return refchr,indexsort,indexname,megam,megarefs,megarefe,refseq,mrnaseq
		
		
		
	def findalltrans(genelist):
			

		t=pd.read_csv('geneinfor.csv', skiprows=0 ,sep ='\t',header=None)
		genename=list(t[1])
		generange=list(t[2])
		trans=list(t[3])
		allgenein={genename[i]:[int(a) for a in generange[i].split(';')] for i in xrange(len(genename)) }		
		genein={genename[i]:trans[i].split('&')[0].split(';') for i in xrange(len(genename))}
			
		list1=pd.read_csv('knownGenefix.txt',sep='\t',header=None)
		transname=list(list1[12].values.flatten())
				
		exons=list(list1[9].values.flatten())
		exone=list(list1[10].values.flatten())
		strand=list(list1[3].values.flatten())
			
		transrecord={transname[i].split('.')[0]:[int(x) for x in sorted(exons[i].split(',')[:-1]+exone[i].split(',')[:-1])] for i in xrange(len(transname))}

		t=pd.read_csv('transinfor.csv', skiprows=0 ,sep ='\t',header=None)
		t=t.iloc[1309574:]
		
		genesname=list(t[2])
		transname=list(t[1])
		genetotrans={genesname[i]:[transname[i]] for i in xrange(len(genesname))}
		allgene={genesname0:[] for genesname0 in genesname}
		genetotrans.update(genein)	
		
		allgene.update(allgenein)

		alltrans={}
		alltransname={}
		for i,genename in enumerate(genelist):
			alltransname[genename]=genetotrans[genename]
			alltrans[genename]=[transrecord[transname] for transname in alltransname[genename]]
		
		return alltransname,alltrans,allgene
	
	print file0		
	t=pd.read_csv(file0, skiprows=1 ,sep =',',header=None)
	
	mstart0=list(t[[8]].values.flatten())
	mend0=list(t[[9]].values.flatten())
	seq0=list(t[[25]].values.flatten())

	mstart0=[[int(a) for a in b.split(';')] for b in mstart0]
        mend0=[[int(a) for a in b.split(';')] for b in mend0]
        for i in xrange(len(seq0)):
                a=len(seq0[i].replace('-',''))
                b=sum(mend0[i])-sum(mstart0[i])
                if a!=b:
					t.drop(t.index[i])




	name0=list(t[[1]].values.flatten())
	ref0=list(t[[2]].values.flatten())
	chrom0=list(t[[3]].values.flatten())
	piecei0=list(t[[12]].values.flatten())
	rstart0=list(t[[13]].values.flatten())
	rend0=list(t[[14]].values.flatten())
	score0=[int(a) for a in list(t[[19]].values.flatten())]
	gstart0=[int(a) for a in list(t[[4]].values.flatten())]
	gend0=[int(a) for a in list(t[[5]].values.flatten())]
	mstart0=list(t[[8]].values.flatten())
	mend0=list(t[[9]].values.flatten())
	seq0=list(t[[25]].values.flatten())
	seq1=list(t[[26]].values.flatten())
	strand0=list(t[[23]].values.flatten())
	strand1=list(t[[24]].values.flatten())
	

	
	tempfolder='temp/{:s}/'.format(outputfile.split('/')[-1])
	try:
		os.mkdir(tempfolder)
	except:
		pass
	
	mstart0=[[int(a) for a in b.split(';')] for b in mstart0]
	mend0=[[int(a) for a in b.split(';')] for b in mend0]
		
	

	refchr,indexsort,genename,megam , megarefs,megarefe,refseq,mrnaseq=combineretro(mstart0,mend0,name0,ref0,chrom0,strand1,piecei0,rstart0,rend0,gstart0,gend0,score0,seq0,seq1)


	alltransname,alltrans,allgenes=findalltrans(list(set(name0)))

	coordinates=[[megam[i]]+alltrans[genename[i]] for i in xrange(len(genename))]
	
	
	allchrname=['chr22.fa', 'chrX.fa', 'chr5.fa', 'chr6.fa', 'chr12.fa', 'chr16.fa', 'chr7.fa', 'chr10.fa', 'chr14.fa', 'chr15.fa', 'chr8.fa', 'chr19.fa', 'chr9.fa', 'chrY.fa', 'chr18.fa', 'chrM.fa', 'chr17.fa', 'chr3.fa', 'chr13.fa', 'chr20.fa', 'chr11.fa', 'chr2.fa', 'chr1.fa', 'chr4.fa', 'chr21.fa']
	
	allchr={}
	for chrname in allchrname:
		
		with open('chroms/'+chrname,mode='r') as f:
			
			allchr[chrname]=''.join(f.read().split('\n')[1:])
			
		f.close()
	
	mrnasize=[[] for a in genename]
	
	for i in xrange(len(mrnasize)):
		
		mrnasize[i]=needleallfiles(i,tempfolder,[megam[i]]+alltrans[genename[i]],allchr[chrom0[indexsort[i][0]]],refseq[i])
		
	del allchr
		
	
	memoryneed=[len(refseq[i])*max([sum([coordinate0[2*x+1]-coordinate0[2*x] for x in xrange(len(coordinate0)/2)])   for coordinate0 in coordinates[i]])  for i in xrange(len(coordinates))]	
			
	m0=16
	record=[0 for a in genename]
	p=mul.Pool(processes=m0)
	for i in xrange(len(record)):
		if memoryneed[i]>90000000:
			print i,' big sequence, skip, size=', memoryneed[i]
			continue
		#findbest(i,tempfolder,[megam[i]]+alltrans[genename[i]],chrom0[indexsort[i][0]],refseq[i],mrnasize[i])			
		record[i]=p.apply_async(findbest,(i,tempfolder,[megam[i]]+alltrans[genename[i]],chrom0[indexsort[i][0]],refseq[i],mrnasize[i]))
		
	p.close()
	p.join()
	scores,bestmrna,highscore=[],[],[]

	unfinished_gene0=range(len(record))
	
	runcycle=0
	k=8
	while len(unfinished_gene0)>0:
		runcycle=runcycle+1
		k=max(1,k-3)	
		p=mul.Pool(processes=k)
		unfinished_gene=[]
		for i in unfinished_gene0:
			
			if runcycle>3 or memoryneed[i]>1300000000:
				print "too many tries, abort",i,genename[i],'\n'
				
				record[i]=[-1,[1.0]+[0.0 for a in alltrans[genename[i]]],[]]

				continue
			try:
				record[i]=record[i].get()
				
			except:			
				unfinished_gene.append(i)
				
				if memoryneed[i]>1300000000.0/k:
					continue

				record[i]=p.apply_async(findbest,(i,tempfolder,[megam[i]]+alltrans[genename[i]],chrom0[indexsort[i][0]],refseq[i],mrnasize[i]))

		
		p.close()
		p.join()
		
		unfinished_gene0=unfinished_gene
		

	allchr={}
	for chrname in allchrname:
		
		with open('chroms/'+chrname,mode='r') as f:
			
			allchr[chrname]=''.join(f.read().split('\n')[1:])
			
		f.close()

	for i,record1 in enumerate(record):
		
		bestmrna.append(record1[0])
		scores.append(record1[1])
		highscore.append(record1[2])
		
		if memoryneed[i]>1300000000.0:
			continue
		
		writecomparefiles(i,tempfolder,[megam[i]]+alltrans[genename[i]],allchr[chrom0[indexsort[i][0]]],refseq[i],record1[0]+1)
		
		
	del allchr
		
	gc.collect()
	
	
	for i,j in enumerate(bestmrna):
		
		
		if j==-1:
		
			exonposix=megam[i]
						
			name1='unknown'
			
			
		else:
				
			exonposix=alltrans[genename[i]][j]
			
			name1=alltransname[genename[i]][j]
		
		
		allposi0=[megam[i]]+alltrans[genename[i]]
		allposi0= '&'.join([';'.join([str(x) for x in a]) for a in allposi0])
				
		out=[[genename[i],name1,ref0[0], chrom0[indexsort[i][0]],piecei0[indexsort[i][0]],megarefs[i],megarefe[i],[allgenes[genename[i]][2*a]  for a in xrange(len(allgenes[genename[i]])/2)],[allgenes[genename[i]][2*a+1] for a in xrange(len(allgenes[genename[i]])/2)],exonposix,strand0[indexsort[i][0]], strand1[indexsort[i][0]],scores[i][j+1],indexsort[i],';'.join(alltransname[genename[i]]),';'.join([str(a) for a in scores[i]]),allposi0]]
	
		output=pd.DataFrame.from_records(out).to_csv(outputfile,sep=',',header=None,mode="a",index=False)

	
		
	return 0



def runcombine(reffiles):
	
	for reffile in reffiles:
		
		file0='result/{:s}_retro.csv'.format(reffile)

		outputfile=file0+'_combine'

		try:
			os.system('mv {:s} {:s}_bak'.format(outputfile,outputfile))
		except:
			pass
	
		combinemrna(file0,outputfile)
	
	
	return 0


	















def findtron(x,y,read,posi,mode):
	
	
	
	def findspan(cordifix,posi0,out=0):

		if posi0==[]:
			if out==0:
				return [],[]
			elif out==1:
				return []
			elif out==2:
				return []

		
		posi=[a-posi0[0] for a in posi0]
			
		#gapcordi=[a.span() for a in re.finditer(r'{:s}+'.format(x),ref)]
				
		#cordifix=sorted(findcodi(query,sorted([a[0] for a in gapcordi]))+[a+1 for a in findcodi(query,sorted([a[1]-1 for a in gapcordi]))])
		
		l0=len(cordifix)
		
		allposi=cordifix+posi
		key0=sorted(range(len(allposi)), key=lambda x: allposi[x])
		
		p0=0
		ps=1
		g0=0
		gs=1
		intron=[]
		exon=[]
		gap0=0
		for i in xrange(len(key0)):
			
			if key0[i]>=l0:
				
				if gs==-1:
					gap0=gap0+posi[key0[i]-l0]-g0	
					g0=posi[key0[i]-l0]			
				if ps==-1:
					exon.append(gap0)
				else:
					intron.append(gap0)
				ps=ps*-1
				gap0=0
				
			else:	
				if gs==1:
					g0=cordifix[key0[i]]
				else:
					
					gap0=gap0+cordifix[key0[i]]-g0
					
				gs=gs*-1
		
		if out==0:
			return exon,intron[1:]
		elif out==1:
			return exon
		elif out==2:
			return intron[1:]
			

	def distractgap(x,query,ref):
		gapcordi=[a.span() for a in re.finditer(r'{:s}+'.format(x),ref)]
		cordifix=[a-query[:a].count('-') for b in gapcordi for a in b]
		return gapcordi,cordifix


	def findstartend(gapcordi,lenth):
		
		if gapcordi==[]:
			return 0,0
			
		
		lastposi=0
		score=0
		eo=-1
		start=0
		for a in gapcordi:
			eo=eo*-1
			if eo>0:
				if a-lastposi>20:
					start=lastposi
					break
				c=4
			else:
				c=1
			score=score+(a-lastposi)*eo*c

			if score<0:
				start=a
				
			lastposi=a

		
			
		gapcordi2=[lenth-a for a in gapcordi][::-1]
		
		
		lastposi=0
		score=0
		eo=-1
		end=0
		for a in gapcordi2:
			eo=eo*-1
			if eo>0:
				if a-lastposi>20:
					end=lastposi
					break
				
				c=4
			else:
				c=1
			score=score+(a-lastposi)*eo*c

			if score<0:
				end=a
				
			lastposi=a


		if start+end>=lenth:
			return -1,-1

		return start,end

	def polish(query,ref,posi0):
		
		l=len(ref)
		
		gapcordi0=[b for a in re.finditer(r'{:s}+'.format('-'),ref) for b in a.span()]
		astart0,aend0=findstartend(gapcordi0, l)
		
		gapcordi1=[b for a in re.finditer(r'{:s}+'.format('-'),query) for b in a.span()]
		astart1,aend1=findstartend(gapcordi1, l)
		
		gapcordi2=[b for a in re.finditer(r'{:s}+'.format('N'),ref) for b in a.span()]
		astart2,aend2=findstartend(gapcordi2, l)
		
		gapcordi3=[b for a in re.finditer(r'{:s}+'.format('N'),query) for b in a.span()]
		astart3,aend3=findstartend(gapcordi3, l)
		astart,aend=max(astart0,astart1,astart2,astart3), max(aend0,aend1,aend2,aend3)
		if -1 in [astart0,astart1,astart2,astart3]:
			return [],-1,-1

		exonsizes=[posi0[2*a+1]-posi0[2*a] for a in xrange(len(posi0)/2)]

		fstream=astart
		for i,exonsize0 in enumerate(exonsizes):
			fstream=fstream-exonsize0
			if fstream<=0:
				break
		if fstream>-10:
			astart=astart-fstream

		fstream=aend
		for i,exonsize0 in enumerate(exonsizes[::-1]):
			fstream=fstream-exonsize0
			if fstream<=0:
				break
		if fstream>-10:
			aend=aend-fstream
			
		
		return gapcordi0,astart,aend
		


	
	if read=='':
		return []
	
		
	alllines=read.split('\n')
	
	#ref: original genome, query:assemblies
	rlines=[a for a in xrange(len(alllines)) if re.match(r'\s*s0',alllines[a])!=None]
	qlines=[a+2 for a in rlines]
	alines=[a+1 for a in rlines]

	if qlines==[]:
		return [],[]
	
	space0=re.search('^\s+',alllines[qlines[0]]).span()[1]
	query=''.join([re.search(r'(?<=[0-9]) ([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
	ref=''.join([re.search(r'(?<=[0-9]) ([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])
	align=''.join([re.search(r'(\s|:)+',alllines[a][space0+3:]).group() for a in alines])[:len(ref)]
	
	
	

	if mode==0:	
		
		lref=len(ref)
		
		#cut both ends of alignments, get better quality
		gapcordi,astart,aend=polish(query,ref,posi[0])
		
		
		if astart==-1:
			print 'no align\n', ref	
			return 0,0,0,0, '','','',[],[]
		
		#find cutted range on mrna
		realrange0=[astart-query[:astart].count('-'),lref-aend-query[:lref-aend].count('-')]
		
		
		qstart0,qend0=realrange0[0], realrange0[1]
		
		#find cutted range on reference geneme
		realrangeabs=cleansame(cutrange(realrange0, [posi[0][2*k] for k in xrange(len(posi[0])/2)],[posi[0][2*k+1] for k in xrange(len(posi[0])/2)]))
		
		if realrangeabs==[]:
			return 0,0,0,0, '','','',[],[]
			
		#find cutted range on assemblly geneme
	
		realassem_range0=[astart-ref[:astart].count('-'),lref-aend-ref[:lref-aend].count('-')]
		
		rstart0,rend0=realassem_range0[0],realassem_range0[1]
		
		#change alignments result

		seq0=query[astart:lref-aend]
		
		seq1=ref[astart:lref-aend]
		
		seq2=align[astart:lref-aend]
			
		
		#findmatch range's cordinates on alignment and mrna
		gapcordi2,cordifix2=distractgap(':',seq0,seq2)
		
		
		#findmatch each exon alignment
		#exon0=findspan(cordifix2,findposi1([a-realrangeabs[0] for a in realrangeabs]),1)
		
		
		#exons=[findcommon(matchrange, posi0) for posi0 in posi]

		return qstart0,qend0,rstart0,rend0, seq0,seq1,seq2,cordifix2,realrangeabs
	
	else:
		print 'findtron', posi	
		gapcordi,cordifix=distractgap('-',query,ref)
				
		intron1=findspan(cordifix, posi,2)
		
		gapcordi,cordifix=distractgap(':',query,align)
			
		exon1=findspan(cordifix, posi,1)

		return exon1,intron1,query,ref,align


def readoutput(file1,file2,file3,posi,mode):
			
	print 'start run\n'	
	cm='stretcher {:s} {:s} -snucleotide2  -gapopen 6 -gapextend 1  {:s}'.format(file1,file2,file3)
	os.system(cm)

	print 'end run\n'
	try:	
		with open(file3, mode='r') as f:
			read=f.read()
		f.close()
	except:

		if mode!=0:
			
			return 0,'',[],[],'','',''
			
		else:
			
			return 0,'',0,0,0,0, '','','',[],[]
			
	try:
		score0=re.search(r'(?<=# Score: )\s*-*[0-9]+\.*[0-9]*', read).group().strip()
	except:
		score0=''
	try:
		inden0=re.search(r'(?<=# Identity: )\s*[0-9]+/[0-9]+', read).group().strip()
	except:
		inden0=''

	if mode!=0:
		
		exon1,intron1,seq3,seq4,seq5=findtron('s0','s1',read,posi,1)

		return score0,inden0,exon1,intron1,seq3,seq4,seq5
		
	else:
		
		qstart0,qend0,rstart0,rend0, seq0,seq1,seq2,matchrange,realrange=findtron('s0','s1',read,posi,0)
		
		return score0,inden0,qstart0,qend0,rstart0,rend0, seq0,seq1,seq2,matchrange,realrange
	

def comparescore(ffile,mfile,ofile,rfile,posi1):
	
	def findsize0(range0):
		
		return sum([range0[2*k+1] for k in xrange(len(range0)/2)])-sum([range0[2*k] for k in xrange(len(range0)/2)])

	def sortpiece(matchrange,range0,originalrange,qstart0,qend0):
		if range0==[] or originalrange==[]:
			return [],[],[]
		
		intronsize=[range0[2*a+2]-range0[2*a+1] for a in xrange(len(range0)/2-1)]
		minusintron=[0]+[a+1 for a in xrange(len(intronsize)) if intronsize[a]<0]+[len(intronsize)+1]
		range0=[range0[2*minusintron[a]:2*minusintron[a+1]] for a in xrange(len(minusintron)-1)]
		
		
		if len(range0)>1:
			
			
			pieceqstart=0
			error=[]
			errorexon=[]
			exoncount=0
			for i,range00 in enumerate(range0):
				
				pieceqend=pieceqstart+sum([range00[2*a+1]-range00[2*a] for a in xrange(len(range00)/2)])
				sizematch=findcommon([pieceqstart,pieceqend], matchrange)
				sizematch=sum([sizematch[2*a+1]-sizematch[2*a] for a in xrange(len(sizematch)/2)])
				
				
				if sizematch<=20:
					
					error.append(i)
					
					errorexon.extend(range(exoncount,exoncount+len(range00)/2))
				
				exoncount=len(range00)/2+exoncount
					
				pieceqstart=pieceqend
				
					
				
		else:
			errorexon=[]	
			error=[]
		
		originalpiece=[[originalrange[2*k],originalrange[2*k+1]] for k in xrange(len(originalrange)/2)]
		
		qposi=0
		exoncount=0
		for a in originalpiece:
			
			qposi=qposi+a[1]-a[0]
			if qstart0>=qposi:
				
				exoncount=exoncount+1
			
			else:
				
				break
		
		qposi=0
		exoncount2=0
		for a in originalpiece:
			
			qposi=qposi+a[1]-a[0]
			
			if qend0>qposi:
				
				exoncount2=exoncount2+1
			
			else:
				
				break
			
		involveexon=range(exoncount,exoncount2+1)
		
		involveexon=[involveexon[a] for a in xrange(len(involveexon)) if a not in errorexon]


		print 'range0',range0		
		return range0,error, involveexon		
	

	
	originalrange=posi1[0]
	
	score1,inden1,qstart0,qend0,rstart0,rend0, seq0,seq1,seq2,matchrange,realrange=readoutput(ofile,mfile,rfile,posi1,0)
	
	if realrange==[]:
		
		return [],[],[],[],0,0,0,0, '','','','','','',[],[]
	
	try:
		inden11=float(inden1.split('/')[0])/float(inden1.split('/')[1])
	except:
		inden11=''
	
	print 'beforesort', realrange	
	orderedpieces,error,involveexon=sortpiece(matchrange,realrange,originalrange,qstart0,qend0)
	
	realrange=[a for i,b in enumerate(orderedpieces) if i not in error for a in b]
	print 'aftersort', realrange
	if realrange==[]:
		return [],[],[],[],0,0,0,0, '','','','','','',[],[]

	posi2=[realrange]+posi1[1:]
		
	score0,inden0,exon1,intron1=[0 for a in orderedpieces],[[0.0,0.0] for a in orderedpieces],[[] for a in orderedpieces],[[] for a in orderedpieces]

	seq3,seq4,seq5='','',''	

	mrnastartpoint=rstart0
	for i,piece0 in enumerate(orderedpieces):

		if i in error:
			mrnastartpoint=mrnastartpoint+sum([piece0[2*a+1]-piece0[2*a]  for a in xrange(len(piece0)/2)])
			continue
		
		startpoint=min(piece0)-min(originalrange)
		
		endpoint=max(piece0)-min(originalrange)
		
		with open(ffile,mode='r') as f:
			read=f.read().split('\n')
			title0=read[0]
			read0=''.join(read[1:])[startpoint:endpoint]
			
				
		f.close()
			
		
		pfile=ffile[:-3]+'_'+str(i)+'.fa'
		
		with open(pfile,mode='w') as f:
			f.write(title0+'\n'+read0)
		f.close()
		
		del read, title0, read0
		
		mrnaendpoint=mrnastartpoint+sum([piece0[2*a+1]-piece0[2*a]  for a in xrange(len(piece0)/2)])
		
		with open(ofile,mode='r') as f:
			read=f.read().split('\n')
			title0=read[0]
			read0=''.join(read[1:])[mrnastartpoint:mrnaendpoint]
						
		f.close()

		mrnastartpoint=mrnaendpoint
		
		qfile=ofile[:-3]+'_'+str(i)+'.fa'
		
		with open(qfile,mode='w') as f:
			f.write(title0+'\n'+read0)
		f.close()
		
		del read, title0, read0
		
		score0[i],inden0[i],exon1[i],intron1[i],seq30,seq40,seq50=readoutput(qfile,pfile,rfile,piece0,1)

		seq3=seq3+seq30
		seq4=seq4+seq40
		seq5=seq5+seq50
	
		try:
			inden0[i]=[float(inden0[i].split('/')[0]),float(inden0[i].split('/')[1])]
		except:
			inden0[i]=[0.0,0.0]
			
		try:
			score0[i]=float(score0[i])
		except:
			score0[i]=0.0
		
		
	
	print realrange,exon1, intron1
	score0,exon1,intron1=sum(score0),[a for b in exon1 for a in b],[a for b in xrange(len(intron1)) for a in intron1[b]+[0] if b not in error][:-1]

	try:
		inden00=sum([a[0] for a in inden0])/sum([a[1] for a in inden0])
	except:
		inden00=0

	inden0='{:f}/{:f}'.format(sum([a[0] for a in inden0]),sum([a[1] for a in inden0]))

	try:
		inden10=float(inden11)/float(inden00)
	except:
		inden10=''	
		
	msize=[findsize0(posi00) for posi00 in posi2]
	
	out=[score0,score1,float(score1)-float(score0),inden0, inden00, inden1,inden11,inden10]
	
	return out,posi2,intron1,exon1,qstart0,qend0,rstart0,rend0, seq0,seq1,seq2.replace(' ','!'),seq3,seq4,seq5.replace(' ','!'),msize,involveexon

	
def fixfile(outputfile,i,gene0,name0, ref0,chrom0,gstart0,gend0,piecei0,rstart0,rend0,strand0,strand1,score1,indexsort,index,scores,posi):
	

	print i,name0,index.split(';')

	if 'unknown' in name0:
		
		posi=[posi[0]]

	else:

		index00=[a.replace('\'','').strip() for a in index.split(';')]
	
		posi=[posi[index00.index(name0.replace('\'','').strip())+1]]

	if posi[0] == []:
		print 'no sequence error\n' 
		return -1

	oriposi=[str(a) for a in posi[0]]
		
	fullsize=max(posi[0])-min(posi[0])
	
	mrnasize=sum([posi[0][2*k+1]-posi[0][2*k] for k in xrange(len(posi[0])/2)])

	if mrnasize==0:
		print 'zero size error\n'
		return -1

	exonnumbers=len(oriposi)/2
	
	tempfolder='temp/{:s}_combine/'.format(outputfile.split('/')[-1].split('_fix')[0])
	
	ffile=tempfolder+'fullgene_{:d}.fa'.format(i)
	
	mfile=tempfolder+'mrna_{:d}.fa'.format(i)
	
	ofile=tempfolder+'refx_{:d}.fa'.format(i)
	
	rfile=tempfolder+'result_{:d}.out'.format(i)
	
	scorefind,posi, intron1,exon1,qstart0,qend0,rstart1,rend1, seq0,seq1,seq2,seq3,seq4,seq5,realsize,involveexon=comparescore(ffile,mfile,ofile,rfile,posi)
	
	if seq0=='':
		
		print "find mistake"+str(i)+'\n'
		return 0
	
	
	exonsize=[posi[0][2*a+1]-posi[0][2*a] for a in xrange(len(posi[0])/2)]
	
	alignmentsize=seq2.count(':')
	
	intronsize=[posi[0][2*a+2]-posi[0][2*a+1] for a in xrange(max(len(posi[0])/2-1,0))]
	
	coverage=realsize[0]*1.000/mrnasize

	if realsize[0]!=0:	

		exonsimi=alignmentsize*1.000/realsize[0]

	else:
		exonsimi=0

	print intron1, intronsize	
	
	try:
		intronlose=sum(intron1)*1.0/sum(intronsize)
	except:
		print 'intron=',intronsize, posi
		intronlose=0.0
	
	eachintronlose=[intron1[a]*1.00/intronsize[a] if intronsize[a]>0 else 0  for a in xrange(len(intron1))]
	
	eachexonsimi=[float(exon1[a])/exonsize[a] if exonsize[a]>0 else 0  for a in xrange(len(exonsize))]
	
	
	try:	
		matchratio=float(scorefind[5].split('/')[0])/float(scorefind[3].split('/')[0])
	except:
		matchratio=-1

	dele=seq0.count('-')
	ins=seq1.count('-')
	genecordi=';'.join(sorted(gstart0+gend0)).replace(' ','')
	if genecordi=='':
		genecordi=';'.join(oriposi)

	if strand1=='+':

		newrstart=rstart0+rstart1
		newrend=rstart0+rend1


	else:
		newrstart=rend0-rstart1
		newrend=rend0-rend1


				
	output= [[i,gene0,name0,ref0,chrom0,genecordi,';'.join(oriposi),';'.join([str(a) for a in posi[0]]),exonnumbers,';'.join([str(a) for a in involveexon]),len(involveexon),float(len(involveexon)*100.0)/exonnumbers,qstart0,qend0,piecei0,newrstart,newrend,fullsize,mrnasize,realsize[0],coverage,score1,seq2.count(':'),seq2.count('!'),seq1.count('-'),seq0.count('-'),strand0,strand1,scorefind[0],scorefind[1],scorefind[2]/mrnasize,scorefind[3],scorefind[4],scorefind[5],scorefind[6],scorefind[7],matchratio,exonsize,exon1,eachexonsimi,exonsimi,intronsize,intron1,eachintronlose,intronlose,indexsort,index,scores,seq0,seq1,seq2,seq3,seq4,seq5]]
	
	lock.acquire()
	output=pd.DataFrame.from_records(output).to_csv(outputfile,sep=',',header=None,mode="a",index=False)
	lock.release()
	
	print "finished evalute"+str(i)+'\n'
	return 0



def evaluate(file0,outputfile):
	
	def readcombine(file0):
		
		t=pd.read_csv(file0, skiprows=0 ,sep =',',header=None)
		
		gene0=list(t[[0]].values.flatten())
		tran0=list(t[[1]].values.flatten())
		ref0=list(t[[2]].values.flatten())
		chrom0=list(t[[3]].values.flatten())
		piecei0=list(t[[4]].values.flatten())
		rstart0=[int(a) for a in list(t[[5]].values.flatten())]
		rend0=[int(a) for a in list(t[[6]].values.flatten())]	
		gstart0=[a[1:-1].split(',') for a in list(t[[7]].values.flatten())]
		gend0=[a[1:-1].split(',') for a in list(t[[8]].values.flatten())]
		allposi=list(t[[16]].values.flatten())
		allposi=[[[int(x) for x in a.split(';')] if a !='' else []  for a in allposi0.split('&')] for allposi0 in allposi]
		score1=[int(a) for a in list(t[[12]].values.flatten())]
		scores=list(t[[15]].values.flatten())
		indexsort=list(t[[13]].values.flatten())
		indexname=list(t[[14]].values.flatten())
		strand0=list(t[[10]].values.flatten())
		strand1=list(t[[11]].values.flatten())
		
		return gene0,tran0,ref0,chrom0,gstart0,gend0,piecei0,rstart0,rend0,strand0,strand1,score1,indexsort,indexname,scores,allposi
	
	gene0,tran0,ref,chrom0,gstart0,gend0,piecei0,rstart0,rend0,strand0,strand1,score1,indexsort,index,scores,allposi= readcombine(file0)
	
	#print tran0,ref,chrom0,gstart0,gend0,piecei0,rstart0,rend0,strand0,strand1,score1,indexsort,index,exonposi
	
	m0=16
	p=mul.Pool(processes=m0)
	
	record0=[0 for a in xrange(len(tran0))]
	print 'size='+str(len(record0))

	for i in xrange(len(record0)):
		#fixfile(outputfile,i,gene0[i],tran0[i],ref[i],chrom0[i],gstart0[i],gend0[i],piecei0[i],rstart0[i],rend0[i],strand0[i],strand1[i],score1[i],indexsort[i],index[i],scores[i],allposi[i])	
		
		record0[i]=p.apply_async(fixfile,(outputfile,i,gene0[i],tran0[i],ref[i],chrom0[i],gstart0[i],gend0[i],piecei0[i],rstart0[i],rend0[i],strand0[i],strand1[i],score1[i],indexsort[i],index[i],scores[i],allposi[i]))
		
	p.close()
	p.join()
	for i,r in enumerate(record0):
		print i
		r=r.get()

	try:
		os.system('rm -r temp')
	except:
		pass

	
	print "finished evalute {:s}\n".format(file0)
	








def runevaluate(reffiles, chromfiles):		

	for reffile in reffiles:
		
		file0='./result/{:s}_retro.csv_combine'.format(reffile)

		outputfile='./result/{:s}_retro.csv_fix'.format(reffile)

		try:
			os.system('mv {:s} {:s}_bak'.format(outputfile,outputfile))
		except:
			pass
		
		with open(outputfile, mode='w') as f:
			header=['sort','Genename', 'Transname_predicted', 'assemblly', 'chr', 'fullgene_codinates', 'mrna_cordinates', 'alignment_cordinates', 'mrna_exons_number', 'involve_exon', 'involve_exon_number', 'percentage_exon_involve', 'qstart_on_mrna', 'qend_on_mrna', 'alignmented_piece_in_assembly', 'assembly_start', 'assembly_end', 'fullsize', 'mrna_size', 'alignment_size', 'coverage', 'multi_score', 'match', 'mismatch', 'deletion', 'insertion', 'strand_on_reference', 'strand_alignment', 'fullgene_score', 'mrna_score', 'score_increase_per_base', 'similarity_fullgene', 'similarity_fullgene_decimal', 'similarity_mrna', 'similarity_mrna_decimal', 'ratio_of_two_similarity', 'ratio_of_match_number', 'exon_sizes', 'each_exon_match','each_exon_similarity' ,'total_exon_similarity', 'Intron_sizes', 'each_intron_deletion', 'each_intron_lose', 'total_intron_lose', 'index_in_last_step', 'all_mrna_of_gene', 'score_of_each_mrna(1st_is_virtul_mrna)', 'trans_seq', 'retro_seq', 'alignment','trans_seq_introned','retro_seq_gapped','trans_introned_align']
			f.write(','.join(header)+'\n')
		f.close()

		evaluate(file0,outputfile)
		os.system('rm -r temp/{:s}'.format(reffile))	


def runbuild(reffiles0):
	for file0 in reffiles0:
		allfiles=os.listdir(wkfolder)
		if file0 not in allfiles:
			cm='gcloud compute copy-files retro-1:/home/michalsakin/data/{:s}  {:s}{:s} --ssh-key-file=retrokey  --zone us-central1-c'.format(file0, wkfolder, file0)
			os.system(cm)
		print wkfolder, file0
		cm='bowtie2-build-l -f --threads 16 {:s}{:s} {:s}refindex/{:s}'.format(wkfolder,file0,savefolder,file0+'_bowtie')
			
		os.system(cm)

def filterfile(f):

	t=pd.read_csv(f,sep=',')

	genecordi=list(t['fullgene_codinates'])

	aligncordi=list(t['alignment_cordinates'])

	exonsimi=list(t['each_exon_similarity'])
	exonmatch=list(t['each_exon_match'])

	intronlose=list(t['each_intron_lose'])
	size=list(t['alignment_size'])
	coverage=list(t['coverage'])
	rstart=list(t['assembly_start'])
	rend=list(t['assembly_end'])

	sp=[]
	filter0=[]
	
	for i in xrange(len(t)):	
		print i

		temp=genecordi[i]
		temp0=re.findall(r'\d+',temp)
		genecordi0=[int(a) for a in temp0]

		if genecordi0==[]:
			continue				
		
		gstart=min(genecordi0)
		gend=max(genecordi0)
		
		try:
			temp=aligncordi[i]
			temp0=re.findall(r'\d+',temp)
		except:
			print i,temp
		aligncordi0=[int(a) for a in temp0]
		astart,aend=min(aligncordi0),max(aligncordi0)
		
		
		if astart<gstart or aend>gend:
			sp.append(i)
			continue
		
		exonsimi0=re.findall(r'\d+\.*\d*',exonsimi[i])

		exonmatch0=re.findall(r'\d+\.*\d*',exonmatch[i])

		if len([j for j in xrange(len(exonsimi0)) if float(exonsimi0[j])>0.9 and float(exonmatch0[j]) >20])<3:
			continue
			
		if len([float(a) for a in re.findall(r'\d+\.*\d*',intronlose[i]) if a>0.8])<2:
			continue
			
		if coverage[i]<0.3:
			continue
			
		if size[i]<100:
			continue
			
		filter0.append(i)

	
	spt=t.iloc[filter0]
	out0=spt.to_csv('filtered/filtered_{:s}'.format(f.split('_fix')[0].split('result/')[1]),sep=',')
	return 0

	


def select(reffiles):

	try:
		os.mkdir('filtered')
	except:
		pass

	for file0 in reffiles:

		print 'filetering'+file0

		resultfile='result/{:s}_retro.csv_fix'.format(file0)
		filterfile(resultfile)



def run():	
	
	genelist=cl.defaultdict(list)
	
	record=[]

	for reffiles0 in reffiles:
		
		try:
			os.mkdir('temp/{:s}'.format(reffiles0))
		except:
			pass
			
		try:
			os.mkdir('result/{:s}'.format(reffiles0))
		except:
			pass			
	if mb==1:
		runbuild(reffiles)
		time.sleep(10)
	if mh==1:
		runhit(reffiles, chromfiles)
		time.sleep(10)

		gc.collect()

		time.sleep(10)
	
	if mt==1:
		transtogene(reffiles,chromfiles)
		gc.collect()
		time.sleep(10)
	
	
	if ma==1:
		runanalyze(reffiles, chromfiles)
		time.sleep(10)
	
	
	if md==1:
		
		runcombine(reffiles)
		gc.collect()
		time.sleep(10)
		
	if me==1:
		runevaluate(reffiles, chromfiles)
		gc.collect()

	if ms==1:
		select(reffiles)
		gc.collect()

	return 0




def outputsummary(genelist):
	
	allretro=0
	t=pd.read_csv('knownGenefix.txt',skiprows=0,sep ='\t',header=None)
	name0=list(t[12])
	del t
	l=len(name0)
	
	findgene=list(set([a for c in genelist for b in c for a in b.keys()]))
	refgenes=[]
	refnums=[]
	for refrecord in genelist:
		
		
		refgenes.append([a for b in refrecord for a in b.keys()])
		
		refnums.append([a for b in refrecord for a in b.values()])
	
	
	allretro=0
	summary=[]
	for findgene0 in findgene:
		summary0=[findgene0]
		
		for i,refrecord in enumerate(genelist):
			
			if findgene0 not in refgenes[i]:
				
				summary0.append(0)
				
				continue
			
			num0=refnums[i][refgenes[i].index(findgene0)]
			
			summary0.append(num0)
			
			allretro=allretro+num0
		
		summary.append(summary0)
	
	out=pd.DataFrame.from_records(summary)
	f = out.to_csv("result/summary.csv",header='genename\t'+'\t'.join(reffiles), mode="w")
	
	return l,len(findgene),allretro


if __name__=='__main__':
	
	
	def findfiles(inputfile,wkfolder):
		
		reffiles=[]
		if inputfile not in ['-1','big' , 'all']:
			reffiles=inputfile.split(';')
			
		if inputfile == 'all':	
			allfiles=os.listdir(wkfolder)
			reffiles=[allfiles[a] for a in xrange(len(allfiles)) if os.path.isfile(wkfolder+allfiles[a]) and re.search('.fasta$',allfiles[a])!=None]
			reffiles0=[allfiles[a] for a in xrange(len(allfiles)) if os.path.isfile(wkfolder+allfiles[a]) and re.search('.fa$',allfiles[a])!=None]
			reffiles=reffiles+reffiles0
		elif inputfile == 'big':  
			
			allfiles=os.listdir(wkfolder)
			reffiles=[allfiles[a] for a in xrange(len(allfiles)) if os.path.isfile(wkfolder+allfiles[a]) and re.search('.fasta$',allfiles[a])!=None]
			reffiles0=[allfiles[a] for a in xrange(len(allfiles)) if os.path.isfile(wkfolder+allfiles[a]) and re.search('.fa$',allfiles[a])!=None]
			reffiles=reffiles+reffiles0
			reffiles=[a for a in reffiles if '_' not in a]

		return reffiles
		
	
	def checkbuild(reffiles):
		
		unbuildfiles=[]

		allbuildfiles=os.listdir('refindex')

		buildfile0=['_'.join(a.split('_')[:-1]) for a in allbuildfiles  if  re.search('.bt2l$',a)!=None]

		for reffile in reffiles:

			if buildfile0.count(reffile)<6:
				
				for end in ['_bowtie.1.bt2l','_bowtie.2.bt2l','_bowtie.3.bt2l','_bowtie.4.bt2l','_bowtie.rev.1.bt2l','_bowtie.rev.2.bt2l']:

					cm='gcloud compute copy-files results:/home/q5476572/refindex/{:s} refindex/{:s}  --zone us-central1-c --ssh-key-file=/home/q5476572/keys/retrokey'.format(reffile+end,reffile+end)
					try:
						a=os.system(cm)
						if a !=0:
							unbuildfiles.append(reffile)
					except:
						unbuildfiles.append(reffile)
		
		unbuildfiles=list(set(unbuildfiles))

		print 'unbuildfiles=',unbuildfiles

		runbuild(unbuildfiles)
		
	
	def checkhit(reffiles):
		
		unhitfiles=[]
		
		finishedhit=os.listdir('result')

		folders=[a for a in finishedhit if os.path.isfile(savefolder+'result/'+a)==False]
		notfinishedhit=[a for a in reffiles if a not in folders]

		for reffile in reffiles:
	

			if reffile in folders:
				try:
					os.system('rm -r result/{:s}_bak'.format(reffile))
				except:
					pass

				os.system('mv result/{:s} result/{:s}_bak'.format(reffile,reffile))
				
			cm='gcloud compute copy-files results:/home/q5476572/hit/{:s} result/{:s}  --zone us-central1-c --ssh-key-file=/home/q5476572/keys/retrokey'.format(reffile,reffile)
			try:
				x=os.system(cm)
				if x!=0:
					unhitfiles.append(reffile)
			except:
				
				unhitfiles.append(reffile)


		for reffile in unhitfiles:

			try:
				os.system('mv result/{:s}_bak result/{:s}'.format(reffile,reffile))
			except:
				pass

				
		unhitfiles=list(set(unhitfiles))

		print unhitfiles
	
		print 'unhitfiles=',unhitfiles		
	
		checkbuild(unhitfiles)
				
		
	
	
	def prepare(reffiles):
		
		exitfiles=os.listdir(wkfolder)

		notinfiles=[a for a in reffiles if a not in exitfiles]

		print 'transport fasta:',notinfiles		

		for file0 in notinfiles:
			
			cm='gcloud compute copy-files retro-1:/home/michalsakin/data/{:s} {:s}  --zone us-central1-c --ssh-key-file=/home/q5476572/keys/retrokey'.format(file0,wkfolder+file0)
			os.system(cm)

		if unfinished==1 or mh==0:
			checkhit(reffiles)


		if   unfinished==0 and   mh==1 and mb==0:
			
			checkbuild(reffiles)
	
		



	print inputfile
	reffiles=findfiles(inputfile,wkfolder)	
	print 'assemblies are:\n'
	print reffiles

	chromfiles=findfiles(inputchrom,chromfolder)

	print 'chroms are:\n'
	
	print chromfiles
	os.chdir(savefolder)
	
	print 'start initiating programs, reading project files\n'
	
	try:
		os.mkdir('result')
	except:
		print "warnning, result folder exist already\n"
		
	try:
		os.mkdir('temp')
	except:
		print "warnning, temp folder exist already\n"
	
	checkmem()	

	prepare(reffiles)

	run()


