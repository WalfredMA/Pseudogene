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

unfinish=0
opts,args=getopt.getopt(sys.argv[1:],"uf:s:")
for op, value in opts:
	if op=='-f':
		inputfile=value
	if op=='-u':
		unfinish=1
	if op=='-s':
		size00=int(value)

allfiles=os.listdir(os.getcwd())
if 'assemblies' in allfiles:	
	assemfolder='/home/q5476572/data/assemblies/'
else:
	assemfolder='/home/michalsakin/data/'

manager = mul.Manager()
seqs = manager.dict()

def cutmrna(seqs,chr0,std, start,end, name,title=0):

	if title==0:
		title=str(chr0)+'_'+str(start)+'_'+str(end)

	if std=='+':
		with open(name, mode='w') as f:
			f.write('>'+title+'\n'+seqs[chr0][start:end])
		f.close()
	else:

		tran=maketrans('ATCGatcg', 'TAGCtagc')

		with open(name,mode='w') as f:
			f.write('>'+title+'\n'+seqs[chr0][end:start][::-1].translate(tran))
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

	query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
	ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])

	truesize=min(len(query)-query.count('-')-query.count('N'),len(ref)-ref.count('-')-ref.count('N'))
	query,ref=polish(query,ref)

	match=len([i for i in xrange(len(query)) if query[i]==ref[i] and query[i] != 'N'])
	mismatch=len([i for i in xrange(len(query)) if query[i]!=ref[i] and query[i] not in ['N','-'] and ref[i] not in ['N','-']])

	gaps=len([b for a in re.finditer(r'{:s}+'.format('-'),query) for b in a.span()])

	gaps2=len([b for a in re.finditer(r'{:s}+'.format('-'),ref) for b in a.span()])

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
		
		alignments_index0=eachgroup_combine(g_index,g_rstart,g_rend,1000)
		
		
		for each_index in alignments_index0:

			
			each_qstart=[qstart[i] for i in each_index]
			
			each_qend=[qend[i] for i in each_index]
			
			
			each_index_q=eachgroup_combine(range(len(each_qstart)),each_qstart,each_qend,200)
			
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






def findreadcoordi(read):
	
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
	
	if read=='':
		return 0,0,0,0
	
	alllines=read.split('\n')
	
	rlines=[a for a in xrange(len(alllines)) if re.match(r'\s*s0',alllines[a])!=None]
	qlines=[a+2 for a in rlines]
	alines=[a+1 for a in rlines]
	
	if qlines==[]:
		return 0,0,0,0
	
	space0=re.search('^\s+',alllines[qlines[0]]).span()[1]
	query=''.join([re.search(r'(?<=[0-9]) ([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
	ref=''.join([re.search(r'(?<=[0-9]) ([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])
	align=''.join([re.search(r'(\s|:)+',alllines[a][space0+3:]).group() for a in alines])[:len(ref)]
		
	
	l=len(ref)
	
	gapcordi0=[b for a in re.finditer(r'{:s}+'.format('-'),ref) for b in a.span()]
	astart0,aend0=findstartend(gapcordi0, l)
	
	gapcordi1=[b for a in re.finditer(r'{:s}+'.format('-'),query) for b in a.span()]
	astart1,aend1=findstartend(gapcordi1, l)
	
	astart,aend=max(astart0,astart1), max(aend0,aend1)
	
	if -1 in [astart0,astart1]:
		return -1,-1,-1,-1
	
	qstart=astart-query[:astart+1].count('-')
	
	qend=aend-query[end:].count('-')
	
	rstart=astart-ref[:astart+1].count('-')
	
	rend=aend-ref[end:].count('-')

	return qstart,qend,rstart,rend

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

	


def lastzcheck(query,ref,out):

	with open(query+'_infor',mode='r') as f:
		infor=f.read().split('_')
	f.close()

	l0=int(infor[0])

	lN=int(infor[1])	
	
	os.system('makeblastdb -in {:s} -dbtype nucl -parse_seqids -out {:s}'.format(ref,ref+'_db'))
	
	
	os.system('blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 1 -max_target_seqs 10'.format(query,ref+'_db',out))

	try:
		t=pd.read_csv(out, sep=',',header=None)
	except:
		t=[]


	t=t.iloc[[i for i in xrange(len(t)) if list(t[4])[i]>(l0-lN)*0.7]]

	return len(t)
	
	



	
def target(file0,f_or_e):
			
	try:
		t=pd.read_csv(file0,sep='\t')
	except:
		return  [],[],[],[],[]

	if len(t)<1:
			
		return [],[],[],[],[]
		
	insertextends=[]
	posis=[]
	chrs=[]
		
	for i in xrange(len(t)):	

		contig=list(t['name2'])[i]

		size0=int(list(t['size1'])[i])


		strand=list(t['strand2'])[i]


		if strand=='-':
			sd=-1
		else:
			sd=1

		if f_or_e=='f':
			if sd==1:
				
				qposi=int(list(t['end1'])[i])
				insertextend=(size0-qposi)
				posi=int(list(t['end2+'])[i])
			else:
				qposi=int(list(t['end1'])[i])
				insertextend=(size0-qposi)
				posi=int(list(t['zstart2+'])[i])
			

		else:
			if sd==1:
				qposi=int(list(t['zstart1'])[i])
				insertextend=qposi
				posi=int(list(t['zstart2+'])[i])
			else:
				qposi=int(list(t['zstart1'])[i])
				insertextend=qposi
				posi=int(list(t['end2+'])[i])

		
		insertextends.append(insertextend)
		posis.append(posi)
		chrs.append(contig)
		
	nmatchs=list(t['nmatch'])
	sizes=list(t['size1'])
		
	return chrs,insertextends,posis,[100*float(nmatchs[i])/int(sizes[i]) for i in xrange(len(nmatchs))],list(t['strand2'])

def findposi(name,svfolder,rangefile,f_or_e):
	
	time=1
	
	while time<10:
		
		os.system('lastz --format=general:name1,zstart1,end1,strand1,size1,name2,zstart2+,end2+,strand2,size2,text1,text2,score,nmatch,nmismatch --notransition   --ambiguous=n --strand=both --ydrop=1000 --match=1,5 --exact=300 --filter=coverage:{:d} --filter=identity:90 {:s} {:s} > {:s}'.format(20/time,svfolder+name+'_'+f_or_e+str(time),rangefile,svfolder+name+'_output'+f_or_e+str(time)))
		
		chrs,extend,posi,scores,strands=target(svfolder+name+'_output'+f_or_e+str(time),f_or_e)
		
		if len(set([a.split('_')[0] for a in chrs]))==1:
			
			break 

		else:
			
			time=time+1
			
			os.system('python denovacut.py -i {:s} -t {:d}  -{:s} '.format(svfolder+name,time,f_or_e))	
	
	return chrs,extend,posi,scores,strands

def findretro(assem,name,svfolder,retrofile,rangefile,size00,seqs):
	
	os.system('blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 1 -max_target_seqs 150'.format(retrofile, rangefile, svfolder+name+'_outputr'))	

	try:
		t=pd.read_csv(svfolder+name+'_outputr',sep=',',header=None)
	except:
		return  []

	if len(t)<1:

		return []



	qstart,qend,rchr,rstart,rend=list(t[6]),list(t[7]),list(t[1]),list(t[8]),list(t[9])

	rstd=['+' if rstart[i]<rend[i] else '-' for i in xrange(len(t))]

	alignments=align_combine(qstart,qend,rchr,rstd,rstart,rend)

	t=pd.DataFrame.from_records(alignments)

	insertextends=[]
	posis=[]
	chrs=[]


	with open(retrofile, mode='r') as f:
		read=''.join(f.read().split('\n')[1:]).strip()
	f.close()

	lr=len(read)
	Ngap=read.count('N')

	del read

	nmatch=list(t[2])

	index=[i for i in xrange(len(t)) if nmatch[i]>(lr-Ngap)*0.8]

	t=t.iloc[index]


	if len(t)<1:

		return []	

	with open(svfolder+name+'_f1', mode='r') as f:
		
		lenf0=''.join(f.read().split('\n')[1:])
		lenf0=len(lenf0)-lenf0.count('N')
	
	f.close()

	with open(svfolder+name+'_e1', mode='r') as f:

		lene0=''.join(f.read().split('\n')[1:])	
		lene0=len(lene0)-lene0.count('N')		

	f.close()


	runends=[1,1]
	if max(lenf0,lene0)<200:

		return 'short'

	results=[]
	for i in xrange(len(t)):


		contig=list(t[0])[i]

		size0=lr

		std=list(t[1])[i]
		
		qstart=min(list(t[3])[i])
		
		qend=max(list(t[4])[i])
		
		score=100*float(list(t[2])[i])/(lr-Ngap)

		strand=list(t[1])[i]

		if strand=='+':
			
			rend0=int(list(t[6])[i])
			
			rend1=rend0+size00/2*3+(lr-qend)
			
			rstart1=int(list(t[5])[i])
			
			rstart0=max(0,int(rstart1)-size00/2*3-qstart)

			#os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format('chroms/'+contig+'.fa', rstart0,rend0,svfolder+name+'_frontcheck'+str(i), contig+'_'+str(rstart0)+'_'+str(rend0)))
			
			#os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format('chroms/'+contig+'.fa', rstart1,rend1,svfolder+name+'_endcheck'+str(i), contig+'_'+str(rstart1)+'_'+str(rend1)))
			
			posi=[contig,'+',rstart1,rend0,qstart,size0-qend]
			
		else:

						
			rstart1=int(list(t[5])[i])

			rstart0=rstart1+size00/2*3+qstart
			
			rend0=int(list(t[6])[i])
			
			rend1=max(0,rend0-size00/2*3-(lr-qend))

			

			#os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s} -r 1'.format('chroms/'+contig+'.fa', rstart0,rend0,svfolder+name+'_frontcheck'+str(i), contig+'_'+str(rstart0)+'_'+str(rend0)))
			
			#os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s} -r 1'.format('chroms/'+contig+'.fa', rstart1,rend1,svfolder+name+'_endcheck'+str(i), contig+'_'+str(rstart1)+'_'+str(rend1)))
			
			posi=[contig,'-',rstart1,rend0,qstart,size0-qend]

		frontfind=-1.0
		endfind=-1.0
	
		if lenf0>200:	

			cutmrna(seqs,contig,strand,rstart0,rend0,svfolder+name+'_frontcheck'+str(i), contig+'_'+str(rstart0)+'_'+str(rend0)) 

			os.system('python cutN.py -f '+svfolder+name+'_f1')

			frontfind=stretcher(svfolder+name+'_f1',svfolder+name+'_frontcheck'+str(i),svfolder+name+'_frontcheckout'+str(i),abs(rstart0-rend0))

		
		if lene0>200:

			cutmrna(seqs,contig,strand,rstart1,rend1,svfolder+name+'_endcheck'+str(i), contig+'_'+str(rstart1)+'_'+str(rend1))

			os.system('python cutN.py -f '+svfolder+name+'_e1')

			endfind=stretcher(svfolder+name+'_e1',svfolder+name+'_endcheck'+str(i),svfolder+name+'_endcheckout'+str(i),abs(rstart1-rend1))
		
		#lr,Ngap,qstartf,qendf,rstartf,rendf,indenff=stretcher(svfolder+name+'_f1',svfolder+name+'_frontcheck'+str(i),svfolder+name+'_frontcheckout'+str(i),1)
		
		#lr,Ngap,qstarte,qende,rstarte,rende,indenee=stretcher(svfolder+name+'_e1',svfolder+name+'_endcheck'+str(i),svfolder+name+'_endcheckout'+str(i),1)

		results.append([contig,strand,rstart1,rend0,qstart,size0-qend,score,frontfind,endfind])
	
		if 1 ==2 :
			if strand=='+':
				sd=1	
			
			else:
				sd=-1
		
			posi[2]=posi[2]-sd*qendf
			
			posi[3]=posi[3]+sd*qstarte
			
			posi[4]=rendf
			
			posi[5]=rstarte
				
			results.append(posi+[indenff,indenee,score,lr,Ngap])
		
	
	return results





def runlastz(arg):


	name,inputfile,chr_ass,s_ass,e_ass,size00=arg[0],arg[1],arg[2],arg[3],arg[4],arg[5]

	global seqs
	

	out=[name,chr_ass,s_ass,e_ass]

	svfolder='temp/{:s}/'.format(inputfile)

	assem=inputfile.split('coordiadd_')[-1].split('_retro.csv')[0]	
			
	retrofile=svfolder+name+'_retro'

	rangefile='db/hg38'

	rangefile2='db_unmask/hg38'
	


	retroaligns=findretro(assem,name,svfolder,retrofile,rangefile,size00,seqs)

	if  retroaligns=='short':
		
		lock1.acquire()

		write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

		write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_short',mode='a',sep=',',header=None,index=False)

		lock1.release()

		return 0

	alteraligns=sorted([a for a in retroaligns if '_' in a[0]],key= lambda x: abs(x[-1]*x[-2]*x[-3]))

	goodalteraligns=[a for a in alteraligns if a[-3]*a[-2]*a[-1]>650000]

	singlealteraligns=[a for a in alteraligns if max(a[-2],a[-1])>90 and min(a[-2],a[-1])<50]

	retroaligns=sorted([a for a in retroaligns if '_' not in a[0]],key= lambda x: abs(x[-1]*x[-2]*x[-3]))

	goodretroaligns=[a for a in retroaligns if a[-3]*a[-2]*a[-1]>650000]

	singleretroaligns=[a for a in retroaligns if max(a[-2],a[-1])>90 and a[-3]*a[-2]*a[-1]<=650000]

	print 'result',name,retroaligns,alteraligns
	if len(retroaligns)>0 and len(goodretroaligns)==1:


		lock1.acquire()

		write=pd.DataFrame.from_records([out+[retroaligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

		lock1.release()

		return 0

	elif len(retroaligns)>0 and len(goodretroaligns)>1:

		lock1.acquire()

		write=pd.DataFrame.from_records([out+[goodretroaligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

		write=pd.DataFrame.from_records([out+['&'.join([';'.join([str(x) for x in a]) for a in goodretroaligns])]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_multi',mode='a',sep=',',header=None,index=False)


		lock1.release()

		return 0
		
	else:
		
		if len(alteraligns)>0 and len(goodalteraligns)==1:


			lock1.acquire()

			write=pd.DataFrame.from_records([out+[alteraligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
			
			write=pd.DataFrame.from_records([out+[alteraligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_alter',mode='a',sep=',',header=None,index=False)

			lock1.release()

		elif len(alteraligns)>0 and len(goodalteraligns)>1:

			lock1.acquire()

			write=pd.DataFrame.from_records([out+[goodalteraligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
			
			write=pd.DataFrame.from_records([out+[goodalteraligns0] for goodalteraligns0 in goodalteraligns]).to_csv('/home/q5476572/data/denova/'+inputfile+'_alter',mode='a',sep=',',header=None,index=False)


			lock1.release()
		



	if len(singleretroaligns+singlealteraligns)>0:

		lock1.acquire()

		write=pd.DataFrame.from_records([out+[(singlealteraligns+singleretroaligns)[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

		write=pd.DataFrame.from_records([out+[singleretroaligns0] for singleretroaligns0 in singleretroaligns+singlealteraligns]).to_csv('/home/q5476572/data/denova/'+inputfile+'_single',mode='a',sep=',',header=None,index=False)


		lock1.release()


	if len(singleretroaligns)==0 and len(retroaligns)>0:

		lock1.acquire()

		write=pd.DataFrame.from_records([out+[retroaligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)


		write=pd.DataFrame.from_records([out+[retroaligns[-1]]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_find',mode='a',sep=',',header=None,index=False)


		lock1.release()

		return 0
		
	elif len(singleretroaligns)==0:
		
		lock1.acquire()

		write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

		write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_notfind',mode='a',sep=',',header=None,index=False)

		lock1.release()

		return 0

	else:
		
		return 0 	

	


	allretrochr=[a[0] for a in retroaligns]

	allretrostd=[a[1] for a in retroaligns]

	allretrostart=[a[2] for a in retroaligns]

	allretroend=[a[3] for a in retroaligns]

	if len(retroaligns)>0:
	
		retroalign=sorted(retroaligns,key=lambda x : x[6]+x[7],reverse=True)[0]

		if min(retroalign[6],retroalign[7])>90:
	
			out.extend(retroalign)
	
			out.append(';'.join(['_'.join([str(x) for x in a]) for a in retroaligns]))

			lock1.acquire()
		
			write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

			lock1.release()

			return 0
	
	
	chrsf,extendf,posif,scoresf,strandf=findposi(name,svfolder,rangefile,'f')
	
	chrse,extende,posie,scorese,strande=findposi(name,svfolder,rangefile,'e')

		
	if min(len(scoresf),len(scorese))<1:
	
		if len(chrsf) >0:
			for i in xrange(len(chrsf)):
				candidates=[j for j in xrange(len(allretrochr)) if allretrochr[j]==chrsf[i] and allretrostd[j]==strandf[i]]
				alldis=[abs(posif[i]-allretrostart[j]) for j in candidates]

				if len(alldis)>0 and min(alldis)<100:
					
					index00=candidates[alldis.index(min(alldis))]

					out.extend(retroaligns[index00])

					break
						
			write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)

			return 0

		if len(chrse) >0:
                        for i in xrange(len(chrse)):
                                candidates=[j for j in xrange(len(allretrochr)) if allretrochr[j]==chrse[i] and allretrostd[j]==strande[i]]
                                alldis=[abs(posie[i]-allretrostart[j]) for j in candidates]

                                if len(alldis)>0 and min(alldis)<100:

                                        index00=candidates[alldis.index(min(alldis))]

                                        out.extend(retroaligns[index00])

                                        break
                                                        
                        write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
			
			return 0

		lock1.acquire()
		
		write=pd.DataFrame.from_records([out+['NA' for i in xrange(12)]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
		
		lock1.release()
		
		return -1
	
	scores=[]
	indexes=[]
	pairs=[]
	for i,chrsf0 in enumerate(chrsf):
		
		matches=[j for j in xrange(len(chrse)) if chrse[j]==chrsf0 and strande[j]==strandf[i]]
		
		if len(matches)>0:
			
			i0=sorted(matches, key=lambda x:abs(posie[x]-posif[i]))[0]
			
			indexes.append(i0)
			
			pairs.append([chrsf[i],strandf[i],posif[i],posie[i0],extendf[i],extende[i0],int(scoresf[i]),int(scorese[i0])])
		
		
	if 	len(pairs)==0:
		
		lock1.acquire()
		
		write=pd.DataFrame.from_records([out+['NA' for i in xrange(12)]]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
		
		lock1.release()
		
		return -1
	
	alldenovascores=[]
	
	for i in xrange(len(pairs)):
		
		chrsf1,strandf1,posif1,posie1,extendf1,extende1=pairs[i][0],pairs[i][1],pairs[i][2],pairs[i][3],pairs[i][4],pairs[i][5]
		
		cutpiecefile='/home/q5476572/data/denova/'+inputfile+'/'+name+'_cutpiece'+str(i)
		
		globalfile='/home/q5476572/data/denova/'+inputfile+'/'+name+'_global'+str(i)

		if chrsf1+'.fa' in os.listdir('chroms'):
			
			if strandf1=='+':
				if posie1-posif1>20:	
					os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format('chroms/'+chrsf1+'.fa', posif1,posie1,cutpiecefile, chrsf1+'_'+str(posif1)+'_'+str(posie1)))
					lr,Ngap,inden00=stretcher(cutpiecefile,retrofile,globalfile)	

				else:
					lr,Ngap,inden00=0,0,0.0

			else:
				if posif1-posie1>20:
				
					os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s} -r 1'.format('chroms/'+chrsf1+'.fa', posif1,posie1,cutpiecefile, chrsf1+'_'+str(posif1)+'_'+str(posie1)))
					lr,Ngap,inden00=stretcher(cutpiecefile,retrofile,globalfile)

				else:

					lr,Ngap,inden00=0,0,0.0
					
		
		alldenovascores.append(inden00)

		pairs[i].extend([inden00,lr,Ngap])
		
		
	highscore=max(alldenovascores)
	
	if highscore<70:
		
		
		longchr_index=[i for i in xrange(len(pairs)) if '_' not in pairs[i][0]]
		
		if longchr_index!=[]:
			
			thebest=sorted(longchr_index, key= lambda x : pairs[x][6]+pairs[x][7])[-1]
		
		else:
			
			thebest=sorted(range(len(pairs)), key= lambda x : pairs[x][6]+pairs[x][7])[-1]

			
		retrofile=svfolder+name+'_retro'
		
		highscore=max(alldenovascores)
		
		chrsf0,chrse0,strand,posif0,posie0,extendf0,extende0=pairs[thebest][0],pairs[thebest][0],pairs[thebest][1],pairs[thebest][2],pairs[thebest][3],pairs[thebest][4],pairs[thebest][5]
				
		conflict_ref,conflictf_ass,conflicte_ass='N','N','N'
		
		if abs(extendf0)>20:
			
			
			os.system('python cutmrna.py -f {:s} -p {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format(assemfolder+assem,chr_ass,int(s_ass)-abs(extendf0),int(s_ass), '/home/q5476572/data/denova/'+inputfile+'/'+name+'_extendf_ass', chr_ass+'_'+str(int(s_ass)-abs(extendf0))+'_'+str(int(s_ass))))
		

			#os.system('python cutmrna.py -f {:s} -s {:d} -e {:d} -o {:s} -t {:s}'.format(chrsf0+'.fa', min(posif0-extendf0,posif0),max(posif0-extendf0,posif0), '/home/q5476572/data/denova/'+inputfile+'/'+name+'_extendf_ref', chrsf0+'_'+str(posif0-extendf0)+'_'+str(posif0)))

		if abs(extende0)>20:

			os.system('python cutmrna.py -f {:s} -p {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format(assemfolder+assem,chr_ass,int(e_ass),int(e_ass+abs(extende0)), '/home/q5476572/data/denova/'+inputfile+'/'+name+'_extende_ass', chr_ass+'_'+str(e_ass)+'_'+str(int(e_ass) + abs(extende0))))

			#os.system('python cutmrna.py -f {:s}  -s {:d} -e {:d} -o {:s} -t {:s}'.format(chrse0+'.fa', min(chrse0+extende0,chrse0),max(chrse0+extende0,chrse0), '/home/q5476572/data/denova/'+inputfile+'/'+name+'_extende_ref', chrse0+'.fa'+'_'+str(chrse0)+'_'+str(chrse0+extende0)))
		
	else:
		
		thebest=alldenovascores.index(highscore)
		
	out.extend(pairs[thebest])
	
	out.append(';'.join(['_'.join([str(x) for x in a]) for a in pairs]))

	os.system('cp {:s} {:s}'.format(retrofile, '/home/q5476572/data/denova/'+inputfile+'/'))
	
	lock1.acquire()
	
	write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_all',mode='a',sep=',',header=None,index=False)
	
	if highscore<70:
				
		write=pd.DataFrame.from_records([out]).to_csv('/home/q5476572/data/denova/'+inputfile+'_denova',mode='a',sep=',',header=None,index=False)
		
	lock1.release()
	
	return 1

try:
        os.system('mkdir temp')
except:
        pass
try:
        os.system('mkdir denova')
except:
        pass

try:
	os.system('mkdir temp/{:s}'.format(inputfile.split('/')[-1]))
except:
	pass
try:
	os.system('mkdir denova/{:s}'.format(inputfile.split('/')[-1]))
except:
	pass


t0=pd.read_csv(inputfile,sep=',')

if unfinish==1 and inputfile.split('/')[-1]+'_all' in os.listdir('/home/q5476572/data/denova/'):

	t1=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_all' ,sep=',',header=None,names=list(range(5)))

	finished=list(t1[0])

	presort=list(t0['sort'])

	preindex=[x for x in xrange(len(presort)) if presort[x] not in finished]

	t0=t0.iloc[preindex]
	
sort=[str(x) for x in list(t0['sort'])]
genename=list(t0['Genename'])


contig=[str(a) for a in list(t0['alignmented_piece_in_assembly'])]
assemble_start=list(t0['assembly_start'])
assemble_end=list(t0['assembly_end'])


assem='_'.join(inputfile.split('/')[-1].split('_')[1:])
if '.fasta' in assem:
	assem=assem.split('.fa')[0]+'.fasta'
else:
	assem=assem.split('.fa')[0]+'.fa'


if unfinish==0:
	os.system('python cutmrnas.py -f {:s} -s {:d}'.format(inputfile, size00))

for file0 in [x for x in  os.listdir('chroms') if '.fa' in x]:

	with open('chroms/'+file0, mode='r') as f:
		read=f.read()
	f.close()
	
	seqs[read.split('\n')[0][1:]]=''.join(read.split('\n')[1:])

	del read

args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]), size00] for i in xrange(len(t0))]
test=0
if test==1:
	map(runlastz,args)

m0=15
p=mul.Pool(processes=m0)

left=p.map_async(runlastz,args)
p.close()
p.join()


left=left.get()



