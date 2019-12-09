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
from string import maketrans

lock1 = mul.Lock()

opts,args=getopt.getopt(sys.argv[1:],"s:f:")
for op, value in opts:
	if op=='-f':
		inputfile=value
	if op=='-s':
		size00=int(value)

allfiles=os.listdir(os.getcwd())
if 'assemblies' in allfiles:	
	assemfolder='/home/q5476572/data/assemblies/'
else:
	assemfolder='/home/michalsakin/data/'

manager = mul.Manager()

seqs = manager.dict()

assemseqs=manager.dict()

def cutmrna(seqs,chr0,std, start,end, name,title=0):
	
	
	if title==0:
		title=str(chr0)+'_'+str(start)+'_'+str(end)


	if std=='+':
		with open(name, mode='w') as f:
			f.write('>'+title+'\n'+seqs[chr0][int(start):int(end)])
		f.close()
	else:

		tran=maketrans('ATCGatcg', 'TAGCtagc')

		with open(name,mode='w') as f:
			f.write('>'+title+'\n'+seqs[chr0][int(end):int(start)][::-1].translate(tran))
		f.close()

	return 0

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

def polish(query,ref):

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
		return '',''

	query=query[astart:len(query)-aend]

	ref=ref[astart:len(ref)-aend]

	return query, ref

def findscore(read):
	if read=='':
		return 0.0


	alllines=[a.strip() for a in read.split('\n')]
	#ref: original genome, query:assemblies
	alllinesindex=[a for a in xrange(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]

	eachline=alllinesindex[::2]
	qlines=[a for a in eachline]
	rlines=[a+2 for a in eachline]
	alines=[a+1 for a in eachline]
	if qlines==[]:
		return 0.0

	try:
		query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
		ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])

	except:
		return 0.0

	truesize=min(len(query)-query.count('-')-query.count('N'),len(ref)-ref.count('-')-ref.count('N'))
	query,ref=polish(query,ref)

	match=len([i for i in xrange(len(query)) if query[i]==ref[i] and query[i] != 'N'])
	mismatch=len([i for i in xrange(len(query)) if query[i]!=ref[i] and query[i] not in ['N','-'] and ref[i] not in ['N','-']])

	gaps=len([b for a in re.finditer(r'{:s}+'.format('-'),query) for b in a.span()])

	gaps2=len([b for a in re.finditer(r'{:s}+'.format('-'),ref) for b in a.span()])

	if truesize==0:
		return 0.0


	return (100.0*match-400*gaps-400*gaps2-80*mismatch)/truesize




def eachgroup_combine(index,rstart,rend,d):

	if len(index)<2:

		return [index]

	sortindex=sorted(range(len(rstart)),key=lambda x: rstart[x])

	index0=[index[i] for i in sortindex]

	rstart0=[rstart[i] for i in sortindex]

	rend0=[rend[i] for i in sortindex]

	e1=rend0[0]

	combine_index=[]
	combine_index0=[index0[0]]
	for i in xrange(1,len(sortindex)):

		dis1=rstart0[i]-e1

		if dis1<d:

			combine_index0.append(index0[i])


		else:

			combine_index.append(combine_index0)

			combine_index0=[index0[i]]

		e1=max(e1,rend0[i])

	combine_index.append(combine_index0)


	return combine_index





def align_combine(qstart,qend,rchr,rstd,rstart,rend):
	
	groups=cl.defaultdict(list)
	
	for i,chr0 in enumerate(rchr):
		
		groups[chr0+':'+rstd[i]].append(i)
				
	alignments=[]
	for key in groups.keys():
		
		g_index=groups[key]
		
		g_rstart=[rstart[i] for i in g_index]
		
		g_rend=[rend[i] for i in g_index]

		std=key[-1]

		if std=='+':

			alignments_index0=eachgroup_combine(g_index,g_rstart,g_rend,10000)

		else:

			alignments_index0=eachgroup_combine(g_index,g_rend,g_rstart,10000)
		
		
		for each_index in alignments_index0:

			
			each_qstart=[qstart[i] for i in each_index]
			
			each_qend=[qend[i] for i in each_index]
			
			
			each_index_q=eachgroup_combine(range(len(each_qstart)),each_qstart,each_qend,20)
			
			each_qstart=[min([each_qstart[x] for x in  index0]) for index0 in each_index_q]
			
			each_qend=[max([each_qend[x] for x in  index0]) for index0 in each_index_q]
			
			chr0=key.split(':')[0]
			
			std0=key.split(':')[1]
			
			if std0=='+':
						
				each_rstart=min([rstart[i] for i in each_index])
				
				each_rend=max([rend[i] for i in each_index])
				
			
			else:
				
				each_rstart=max([rstart[i] for i in each_index])
				
				each_rend=min([rend[i] for i in each_index])
				
			size0=sum([abs(each_qend[i]-each_qstart[i]) for i in xrange(len(each_qend))])
			
			alignments.append([chr0,std0,size0,each_qstart,each_qend,each_rstart,each_rend])
	
	return alignments

def stretcher(query,ref,globalfile,size0):

	os.system('python cutN.py -f '+query)

	with open(query+'_infor',mode='r') as f:
		infor=f.read().split('_')
	f.close()

	lr=int(infor[0])

	Ngap=int(infor[1])


	cm='stretcher {:s} {:s} -snucleotide2  -gapopen 8  -gapextend 1  {:s}'.format(query,ref,globalfile)

	try:
		os.system(cm)

	except:
		pass

	try:
		with open(globalfile, mode='r') as f:
			read=f.read()
		f.close()

	except:

		read=''

	score=findscore(read)


	return max(-2.0,min(100.0,score))




def blastncheck(query,ref):
	
	os.system('python cutN.py -f '+query)

	with open(query+'_infor',mode='r') as f:
		infor=f.read().split('_')
	f.close()

	lr=int(infor[0])

	Ngap=int(infor[1])
	
	os.system('blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 1 -max_target_seqs 150'.format(query, ref, query+'out'))	

	try:
		t=pd.read_csv(query+'out',sep=',',header=None)
	except:
		return  [],0

	if len(t)<1:

		return [],0



	qstart,qend,rchr,rstart,rend=list(t[6]),list(t[7]),list(t[1]),list(t[8]),list(t[9])

	rstd=['+' if rstart[i]<rend[i] else '-' for i in xrange(len(t))]

	alignments=align_combine(qstart,qend,rchr,rstd,rstart,rend)

	t=pd.DataFrame.from_records(alignments)

	nmatch=list(t[2])

	index=[i for i in xrange(len(t)) if nmatch[i]>(lr-Ngap)*0.5]

	t=t.iloc[index]


	if len(t)<1:

		return [],0	


	results=[]
	for i in xrange(len(t)):


		contig=list(t[0])[i]

		size0=lr

		std=list(t[1])[i]
		
		qstart=min(list(t[3])[i])
		
		qend=max(list(t[4])[i])
		
		score=100*float(list(t[2])[i])/(lr-Ngap)

		strand=list(t[1])[i]
		
		rstart1=int(list(t[5])[i])
		
		rend1=int(list(t[6])[i])

				
		results.append([contig,strand,rstart1,rend1,qstart,size0-qend,score])
	
		
	
	return results,Ngap






def runblastn(arg):
	
	global seqs,assemseqs

	name,inputfile,contig,astart,aend,rchr, rstd, rstart,rend,badend,badscore,size00,table=arg[0],arg[1],arg[2],arg[3],arg[4],arg[5],arg[6],arg[7],arg[8],arg[9],arg[10],arg[11],arg[12]

	if rstd=='+':
		sd=1
	else:
		sd=-1

	svfolder='temp/{:s}/'.format(inputfile)

	rangefile='db/hg38'
	
	if max(astart,len(assemseqs[contig])-aend)<200:
		
		lock1.acquire()
		
		write=table.to_csv('/home/q5476572/data/denova/'+inputfile+'_short',mode='a',sep=',',header=None,index=False)
		
		lock1.release()
		
		return ''

	f1m,b1m='a','a'
	finish=2
	while finish>0:

		if badend==1:

			with open(svfolder+name+'f1m', mode='w') as f:
			
				f.write('>'+str(contig)+'_'+str(max(0,astart-300000-size00))+'_'+str(max(0,astart-300000-size00)+size00)+'\n'+assemseqs[contig][max(0,astart-300000-size00):max(0,astart-300000-size00)+size00])
			
			f.close()

			cutmrna(seqs,rchr, rstd, max(0,rstart-sd*400000),rend, svfolder+name+'f1mref')

			f1m=stretcher(svfolder+name+'f1m',svfolder+name+'f1mref',svfolder+name+'f1mout',0)
			
			if finish==2:
				newscore=f1m	

		if badend==-1:
			
			with open(svfolder+name+'e1m', mode='w') as f:
			
				f.write('>'+str(contig)+'_'+str(min(aend+300000+size00,len(assemseqs[contig]))-size00)+'_'+str(min(aend+300000+size00,len(assemseqs[contig])))+'\n'+assemseqs[contig][min(aend+300000+size00,len(assemseqs[contig]))-size00:min(aend+300000+size00,len(assemseqs[contig]))])
				
			f.close()


			cutmrna(seqs,rchr, rstd, rstart,rstart+sd*400000, svfolder+name+'e1mref')

			b1m=stretcher(svfolder+name+'e1m',svfolder+name+'e1mref',svfolder+name+'e1mout',0)


			if finish==2:
				newscore=b1m

		if newscore>50:

			finish=0

		else:
			finish=finish-1

			badend=badend*-1


	print name,newscore


	if newscore<50:

		table.iloc[0,4]=table.iloc[0,4]+'@[{:f},{:f}]'.format(float(f1m), float(b1m))

		lock1.acquire()
		
		write=table.to_csv('/home/q5476572/data/denova/'+inputfile+'_singlenega',mode='a',sep=',',header=None,index=False)
		
		lock1.release()

		return ''

	else:

		table.iloc[0,4]=table.iloc[0,4]+'@[{:f}]'.format(float(newscore))

		if '_' in rchr:

			lock1.acquire()

			write=table.to_csv('/home/q5476572/data/denova/'+inputfile+'_singlealterposi',mode='a',sep=',',header=None,index=False)
		
			lock1.release()

			return ''

		else:
	
			lock1.acquire()

			write=table.to_csv('/home/q5476572/data/denova/'+inputfile+'_singleposi',mode='a',sep=',',header=None,index=False)

			lock1.release()
		
			return name








try:
        os.system('mkdir temp')
except:
        pass
try:
        os.system('mkdir denova2')
except:
        pass

try:
	os.system('mkdir temp/{:s}'.format(inputfile.split('/')[-1]))
except:
	pass
try:
	os.system('mkdir denova2/{:s}'.format(inputfile.split('/')[-1]))
except:
	pass

try:
	os.system('rm denova/{:s}'.format(inputfile.split('/')[-1]+'_singleposi'))

except:

	pass

try:
	os.system('rm denova/{:s}'.format(inputfile.split('/')[-1]+'_singlenega'))

except:

	pass

file0=inputfile.split('/')[-1]

t1=pd.read_csv('denova/'+file0+'_single',sep=',',header=None)

sort1=list(t1[0])

infor=[a[1:-1].split(',') for a in list(t1[4])]


front=[float(a[-2].strip()) for a in infor]
back=[float(a[-1].strip()) for a in infor]


single=[i for i in xrange(len(front)) if max(front[i],back[i])>90 and  min(front[i],back[i])!=-1.0]

t1=t1.iloc[single]
t1=t1.sort_values(by=[0])
infor=[a[1:-1].split(',') for a in list(t1[4])]


sort=[str(x) for x in list(t1[0])]
chrom=[a[0].strip()[1:-1] for a in infor]
std=[a[1].strip()[1:-1] for a in infor]
start=[int(a[2].strip()) for a in infor]
end=[int(a[3].strip()) for a in infor]

badend=[1 if back[i]>front[i] else -1 for i in xrange(len(back))]
badscore=[min(back[i],front[i]) for i in xrange(len(back))]


t0=pd.read_csv(inputfile,sep=',')
presort=[str(x) for x in list(t0['sort'])]
preindex=[presort.index(x)  for x in sort]

t0=t0.iloc[preindex]
genename=list(t0['Genename'])


contig=[str(a) for a in list(t0['alignmented_piece_in_assembly'])]
assemble_start=list(t0['assembly_start'])
assemble_end=list(t0['assembly_end'])




assem='_'.join(file0.split('_')[1:])
if '.fasta' in assem:
	assem=assem.split('.fa')[0]+'.fasta'
else:
	assem=assem.split('.fa')[0]+'.fa'


if assem in os.listdir(assemfolder):
	
	with open(assemfolder+assem, mode='r') as f:
		read=f.read().split('>')[1:]
	f.close()
	for read0 in read:
		
		assemseqs[read0.split('\n')[0].split(' ')[0]]=''.join(read0.split('\n')[1:])

	del read


#os.system('python cutmrnas.py -f {:s} -s {:d}'.format(inputfile, 30000))




for file0 in [x for x in  os.listdir('chroms') if '.fa' in x]:

	with open('chroms/'+file0, mode='r') as f:
		read=f.read()
	f.close()
	
	seqs[read.split('\n')[0][1:].split(' ')[0]]=''.join(read.split('\n')[1:])

	del read
test=0
if test==1:
	args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]),chrom[i],std[i],start[i],end[i],badend[i],badscore[i],size00,t1.iloc[[i]]] for i in xrange(len(t1)) if '_'  in chrom[i]]

	left=map(runblastn,args)
	exit()

m0=15
p=mul.Pool(processes=m0)

args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]),chrom[i],std[i],start[i],end[i],badend[i],badscore[i],size00,t1.iloc[[i]]] for i in xrange(len(t1)) if '_' not  in chrom[i]]
left=p.map_async(runblastn,args)
p.close()
p.join()


left=left.get()
left=list(set(left))

m0=15
p=mul.Pool(processes=m0)

args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]),chrom[i],std[i],start[i],end[i],badend[i],badscore[i],size00,t1.iloc[[i]]] for i in xrange(len(t1)) if '_' in chrom[i] and str(sort[i]) not in left]
left=p.map_async(runblastn,args)
p.close()
p.join()



left=left.get()





