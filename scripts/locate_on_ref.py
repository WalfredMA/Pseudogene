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

opts,args=getopt.getopt(sys.argv[1:],"c:f:")
for op, value in opts:
	if op=='-f':
		inputfile=value

allfiles=os.listdir(os.getcwd())
if 'assemblies' in allfiles:	
	assemfolder='/home/q5476572/data/assemblies/'
else:
	assemfolder='/home/michalsakin/data/'

manager = mul.Manager()
seqs = manager.dict()
assemseqs = manager.dict()

def distinguish(x,y):

	if x==y and x!='N':

		return 1
		
	elif x=='N' or y=='N':
		
		return 2

	elif y=='-':

		return 3


	elif x=='-':
	
		return 0

	else:

		return -1



def findgap(querygap):
	
	if len(querygap)<2:
		
		return [[x[0],x[1]] for x in querygap]
	
	gapsize=[x[1]-x[0] for x in querygap]
	
	seqsize=[querygap[i+1][0]-querygap[i][1] for i in xrange(len(querygap)-1)]
	
	
	insertion=[]
	s0=querygap[0][0]
	e0=querygap[0][1]
	for i,gapsize0 in enumerate(gapsize[:-1]):
		
		if seqsize[i]>min(gapsize0+gapsize[i+1],50):
			
			insertion.append([s0,e0])
			s0=querygap[i+1][0]
			e0=querygap[i+1][1]
			
		else:
			
			e0=querygap[i+1][1]
			
	insertion.append([s0,e0])

	return insertion

def cutmrna(seqs0,chr0,std, start,end, name,title=0):
	if title==0:

		title=str(chr0)+'_'+str(start)+'_'+str(end)


	if std=='+':
		with open(name, mode='w') as f:
			f.write('>'+title+'\n'+seqs0[chr0][start:end])
		f.close()

	else:
		tran=maketrans('ATCGatcg', 'TAGCtagc')
		with open(name,mode='w') as f:
			f.write('>'+title+'\n'+seqs0[chr0][end:start][::-1].translate(tran))
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
			if a-lastposi>100:
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
			if a-lastposi>100:
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
	
	
	default=[[-1,-1,-1,-1,-1,-1,-1],0.0,0.0,0]
	
	if read=='':
		return default


	alllines=[a.strip() for a in read.split('\n')]
	#ref: original genome, query:assemblies
	alllinesindex=[a for a in xrange(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]

	eachline=alllinesindex[::2]
	qlines=[a for a in eachline]
	rlines=[a+2 for a in eachline]
	alines=[a+1 for a in eachline]
	if qlines==[]:
		return default

	query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
	ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])

	true_qsize=len(query)-query.count('-')-query.count('N')
	
	true_rsize=len(ref)-ref.count('-')-ref.count('N')

	truesize=min(true_qsize,true_rsize)
	
	if truesize==0:
		
		return default
	
	query,ref=polish(query,ref)
	
	if len(query)==0:
		
		return default
		
	
	scores=[distinguish(query[i],ref[i]) for i in xrange(len(query))]
	
	match=scores.count(1)
	
	mismatch=scores.count(-1)
	
	Ngap=scores.count(2)
	
	deletion=scores.count(0)
	
	insertion=scores.count(3)
	
	align_size=len(query)-Ngap
	
	if align_size<50:
		
		return default
	
	
	dele_count=len([b for a in re.finditer(r'{:s}+'.format('-'),query) for b in a.span()])

	insert_count=len([b for a in re.finditer(r'{:s}+'.format('-'),ref) for b in a.span()])
	
	scores=[100*match/align_size,100*mismatch/align_size,100*Ngap/len(query),100*deletion/align_size,100*insertion/align_size,dele_count,insert_count]
	
	global_score=100.0*match-320*dele_count-320*insert_count-100*mismatch
	
	local_score=100.0*match-320*dele_count-320*insert_count-100*mismatch-100*insertion
	

	return [scores,global_score,local_score,align_size]







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





def align_combine(qstart,qend,rchr,rstd,rstart,rend,size0=10000):
	
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

			alignments_index0=eachgroup_combine(g_index,g_rstart,g_rend,size0)

		else:
			
			alignments_index0=eachgroup_combine(g_index,g_rend,g_rstart,size0)	

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

def stretcher(query,ref,globalfile,mode='global'):
	
	
	with open(query,mode='r') as f:
		read0=''.join(f.read().split('\n')[1:])
	f.close()

	lenq=len(read0)-read0.count('N')

	del read0

	with open(ref,mode='r') as f:
		read0=''.join(f.read().split('\n')[1:])
	f.close()

	lenr=len(read0)-read0.count('N')

	del read0
	
	
	if min(lenq,lenr)<100:
		
		if mode=='local':
		
			return [-1,-1,-1,-1,-1,-1,-1],0.0,0,lenq, lenr
		
		else:
			
			return 0.0
	

	cm='stretcher {:s} {:s} -snucleotide2  -gapopen 32  -gapextend 0  {:s}'.format(query,ref,globalfile)

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
		
	
	if read=='':
		
		if mode=='local':
		
			return [-1,-1,-1,-1,-1,-1,-1],0.0,0,lenq, lenr
		
		else:
			
			return 0.0


	if mode=='score':
		
		try:
			score=int(re.search(r'(?<=# Score: )\s*-*[0-9]+\.*[0-9]*', read).group().strip())
		except:
			score=0.0

		return score
	
	

	align_infor=findscore(read)
		
	infor=align_infor[0]

	globalscore=align_infor[1]
	
	localscore=align_infor[2]
	
	align_size=align_infor[3]
	

	if mode=='local':
		

		return infor, localscore, align_size, lenq, lenr
		
		
	elif mode=='global':


		return globalscore/max(1000,min(lenq, lenr))


def examposi(ref0,query0,std,rseq,qseq,name0,fsize=-1,esize=-1):


	cutmrna(query0,qseq[0],'+',max(0,qseq[1]-100000),qseq[1],name0+'qf')
	
	cutmrna(query0,qseq[0],'+',qseq[2],qseq[2]+100000,name0+'qe')
	
	if std=='+':

		cutmrna(ref0,rseq[0],'+',max(0,rseq[1]-100000),rseq[1], name0+'rf')
		cutmrna(ref0,rseq[0],'+',rseq[2],rseq[2]+100000, name0+'re')

	else:

		cutmrna(ref0,rseq[0],'-',rseq[1]+100000,rseq[1], name0+'rf')

		cutmrna(ref0,rseq[0],'-',rseq[2],max(0,rseq[2]-100000), name0+'re')

	if fsize>=0 and fsize<100:

		score1=0.0

	else:

		score1=stretcher(name0+'qf',name0+'rf',name0+'of','global')
		

	if esize>=0 and esize<100:

		score2=0.0

	else:

		score2=stretcher(name0+'qe',name0+'re',name0+'oe','global')


	return max(0.0,score1)+max(0.0,score2)



	
def findcodi(seq,posi):
	if posi==[]:
		return []

	
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

				
def findglobalalign(read,s,e):
	
	if read=='':
		
		return [-1,-1,-1,-1,0.0,-1,-1,-1,-1,-1]


	
	alllines=[a.strip() for a in read.split('\n')]
	#ref: original genome, query:assemblies
	alllinesindex=[a for a in xrange(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]

	eachline=alllinesindex[::2]
	qlines=[a for a in eachline]
	rlines=[a+2 for a in eachline]
	alines=[a+1 for a in eachline]
	if qlines==[]:
		return [-1,-1,-1,-1,0.0,-1,-1,-1,-1,-1]

	query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines])
	ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines])

	querygap=[a.span() for a in re.finditer(r'-+',query)]
	querygap=findgap(querygap)	
	gapinfor=[]
	querycodi=0
	
	lastgap=0
	for x in querygap:
	
		querycodi=querycodi+x[0]-lastgap
		
		refcodi=[x[0]-ref[:x[0]].count('-'),x[1]-ref[:x[1]].count('-')]
		
		gapinfor.append([querycodi,refcodi[0],refcodi[1],x[0],x[1]])
		
		lastgap=x[1]
	
	refcodi=findcodi(ref,[s,e])	
	
	denovogap=[x for x in gapinfor if max(x[2],x[1],s,e)-min(x[2],x[1],s,e)-abs(e-s)-abs(x[2]-x[1])<0]
	
	aligngaps=min([x[3] for x in denovogap]+[refcodi[0]])

	aligngape=max([x[4] for x in denovogap]+[refcodi[1]])
	
	frontqseq=query[:aligngaps]
	
	endqseq=query[aligngape:]
	
	gapcordi0=[b for a in re.finditer(r'{:s}+'.format('-'),frontqseq) for b in a.span()]
	astart0,aend0=findstartend(gapcordi0,aligngaps)


	if -1 in [astart0, aend0]:
		
		aend0=0
	
	realaligngaps=aligngaps-aend0
	
	gapcordi0=[b for a in re.finditer(r'{:s}+'.format('-'),endqseq) for b in a.span()]
	astart0,aend0=findstartend(gapcordi0, len(ref)-aligngape)
	
	if -1 in [astart0, aend0]:

		astart0=0


	realaligngape=aligngape+astart0
	

	querygaps=realaligngaps-query[:realaligngaps].count('-')
	
	querygape=realaligngape-query[:realaligngape].count('-')
	

	refgaps=realaligngaps-ref[:realaligngaps].count('-')

	refgape=realaligngape-ref[:realaligngape].count('-')

	querygapseq=query[refcodi[0]:refcodi[1]]
	
	refgapseq=ref[refcodi[0]:refcodi[1]]
	


	gaps=len([b for a in re.finditer(r'{:s}+'.format('-'),querygapseq) for b in a.span()])

	gaps2=len([b for a in re.finditer(r'{:s}+'.format('-'),refgapseq) for b in a.span()])


	gapmatch=[distinguish(querygapseq[x],refgapseq[x]) for x in xrange(len(querygapseq))]

	if len(gapmatch)-gapmatch.count(0)!=0:

		quality=max(0.0,(100.0*(gapmatch.count(1)+gapmatch.count(2))-320*gaps-320*gaps2-100*gapmatch.count(-1)))/(len(gapmatch)-gapmatch.count(0))

	else:

		quality=100.0

	noinsertionsize=len(gapmatch)-gapmatch.count(2)-gapmatch.count(3)

	if noinsertionsize!=0:

		insertion=100.0*gapmatch.count(0)/noinsertionsize
		
		deletion=100.0*gapmatch.count(3)/noinsertionsize

		mismatch=100.0*gapmatch.count(-1)/noinsertionsize

		match=100.0*gapmatch.count(1)/noinsertionsize

		Nsize=100.0*gapmatch.count(2)/(len(gapmatch)-gapmatch.count(3))

	else:

		insertion=0.0

		mismatch=0.0

		deletion=0.0

		match=100.0
	
		Nsize=100.0

	realevent=[querygaps, querygape, refgaps,  refgape,quality,insertion,deletion,mismatch,Nsize, match]	
	
	return realevent


		
	


def findinsertion(svfolder,inputfile,name,index,ref,lf,lfr,rn=0):
	
	score1=stretcher(ref,'denova2/{:s}/'.format(inputfile)+name+'_fullseq','denova2/{:s}/'.format(inputfile)+name+'_fullseqout'+str(index),'score')
	
	score2=stretcher(ref,'denova2/{:s}/'.format(inputfile)+name+'_fe','denova2/{:s}/'.format(inputfile)+name+'_feout','score')
	
	try:
		with open('denova2/{:s}/'.format(inputfile)+name+'_fullseqout'+str(index), mode='r') as f:
			read=f.read()
		f.close()

	except:

		read=''

	realevent=findglobalalign(read,lf,lfr)

	return score2-score1-5*rn,realevent



def findretro(svfolder,ref,name,i,chr0,std,cutstart,cutend,seqs):


	retrofind,retrosize,retroN=blastncheck(svfolder+name+'_retro',ref+'db',ref+'out'+str(i),50,50)


	if len(retrofind)==0:

		return [std,cutstart,cutend,0.0,0.0,'-',0.0,0.0,100*retroN/retrosize,100.0,0.0]

	highscore=-10000000000000

	for j,retrofind0 in enumerate(retrofind):
		
		if retrofind0[1]=='-':
			
			frame='shift'
		
		else:
			
			frame='-'

		if (std=='+' and retrofind0[1]=='+') or (std=='-' and retrofind0[1]=='-'):

			retrofind0[1]='+'

		else:

			retrofind0[1]='-'


		if std=='+':

			retrostart=cutstart+retrofind0[2]

			retroend=cutstart+retrofind0[3]

		else:
			retrostart=cutstart-retrofind0[2]

			retroend=cutstart-retrofind0[3]


		
		cutmrna(seqs,chr0,retrofind0[1],retrostart,retroend, ref+'cut'+str(i)+'_'+str(j))


		infor,score,alignsize,qsize,rsize=stretcher(ref+'cut'+str(i)+'_'+str(j),svfolder+name+'_retro',ref+'out'+str(i)+'_'+str(j),'local')

		
		if score>highscore:

			highscore=score
			
			out=[retrofind0[1],retrostart,retroend,100.0*alignsize/(retrosize-retroN),highscore/max(1,alignsize),frame]+infor[:5]
			
	return out




def blastncheck(query,ref,out,size00=500,csize=500):
	
	with open(query,mode='r') as f:
		infor=''.join(f.read().split('\n')[1:])
	f.close()

	lr=len(infor)

	Ngap=infor.count('N')

	if lr-Ngap<=csize:

		return  [],lr,Ngap

	
	if size00<500:	
		a=0
	
	os.system('blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 4 -max_target_seqs 150'.format(query, ref, out))	

	try:
		t=pd.read_csv(out,sep=',',header=None)
	except:
		return  [],lr,Ngap

	if len(t)<1:

		return [],lr,Ngap
		
	
	prestart,preend=list(t[8]),list(t[9])
	
	t=t.iloc[[i for i in xrange(len(prestart)) if abs(preend[i]-prestart[i])>size00]]

	if len(t)<1:

		return [],lr,Ngap


	qstart,qend,rchr,rstart,rend=list(t[6]),list(t[7]),list(t[1]),list(t[8]),list(t[9])

	rstd=['+' if rstart[i]<rend[i] else '-' for i in xrange(len(t))]

	alignments=align_combine(qstart,qend,rchr,rstd,rstart,rend,size00*200)

	t=pd.DataFrame.from_records(alignments)

	nmatch=list(t[2])

	index=[i for i in xrange(len(t)) if nmatch[i]>max((lr-Ngap)*0.2,csize)]

	t=t.iloc[index]

	if len(t)<1:

		return [],lr,Ngap	


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
	
		
	return results,lr,Ngap


def runblastn(arg):
	
	global seqs,assemseqs

	name,inputfile,contig,start,end=arg[0],arg[1],arg[2],arg[3],arg[4]

	print 'running', name
	
	svfolder='temp/{:s}/'.format(inputfile)
	
	rangefile='db/hg38'
	
	with open(svfolder+name+'_retro',mode='r') as f:

		read0=f.read()
		retroseq=''.join(read0.split('\n')[1:])

	f.close()

	with open('denova2/{:s}/'.format(inputfile)+name+'_retro',mode='w') as f:

		f.write(read0)

	f.close()

	del read0


	with open(svfolder+name+'_f1',mode='r') as f:

		read0=f.read()
		fronttitle=read0.split('\n')[0]
		frontseq=''.join(read0.split('\n')[1:])

	f.close()

	del read0

	with open(svfolder+name+'_e1',mode='r') as f:

		read0=f.read()
		endtitle=read0.split('\n')[0]
		endseq=''.join(read0.split('\n')[1:])

	f.close()

	del read0

	codis=[int(fronttitle.split('_')[-2]),int(fronttitle.split('_')[-1]),int(endtitle.split('_')[-2]), int(endtitle.split('_')[-1])]

	contig0='_'.join(fronttitle.split('_')[:-2])[1:]
	fulls=codis[0]
	fulle=codis[-1]

	with open('denova2/{:s}/'.format(inputfile)+name+'_fullseq',mode='w') as f:

		f.write('>{:s}_{:d}_{:d}\n'.format(contig0,fulls,fulle)+frontseq+retroseq+endseq)

	f.close()

	with open('denova2/{:s}/'.format(inputfile)+name+'_fe',mode='w') as f:

		f.write('>{:s}_{:d}_{:d}\n'.format(contig0,fulls,fulle)+frontseq+endseq)

	f.close()

	lf,lr,le=len(frontseq),len(retroseq),len(endseq)
	fn,rn,en=frontseq.count('N'), retroseq.count('N') ,endseq.count('N')
	lfr=lf+lr
	
	if min(lf-fn,le-en)<100:
	
		lock1.acquire()
		
		write=pd.DataFrame.from_records([[name,0.0,lr,lf,le,rn,fn,en]]).to_csv('/home/q5476572/data/denova2/'+inputfile+'_unknown',mode='a',sep='\t',header=None,index=False)
		
		lock1.release()
	
		return 0 

	
	del frontseq,retroseq,endseq

	frontfind,lf0,frontN=blastncheck(svfolder+name+'_f1',rangefile,svfolder+name+'_f1out')
	
	endfind,le0,endN=blastncheck(svfolder+name+'_e1',rangefile,svfolder+name+'_e1out')
	
	frontfind=sorted(frontfind,key=lambda x :x[-1],reverse=True)

	with open(svfolder+name+'_retro',mode='r') as f:
		retroseq=''.join(f.read().split('\n'))
		retroN=retroseq.count('N')
		retrosize=len(retroseq)
	f.close()
	

	allout=[]

	posiscore=[]
	negascore=[]

	allscore=[]
	mainscore=[]
	allsamechre=[]
	allanchorscore=[]
	find=0

	for i,frontfind0 in enumerate(frontfind):

		chr0,std,s0,e0=frontfind0[0],frontfind0[1],frontfind0[2],frontfind0[3]
		
		samechre=[a for a in endfind if a[0]==chr0]

		dis=[]
		for samechre0 in samechre:

			cutstart=min(s0,e0,samechre0[2],samechre0[3])

			cutend=max(s0,e0, samechre0[2],samechre0[3])
				
			dis0=cutend-cutstart-abs(s0-e0)-abs(samechre0[2]-samechre0[3])

			dis.append(dis0)
				
		samechreindex=sorted([x for x in xrange(len(samechre)) if dis[x]<100000],key=lambda x : dis[x])

		samechre=[samechre[x] for x in samechreindex]	

		if len(samechre)==0:

			continue

		bothscore=frontfind0[6]+samechre[0][6]

		allsamechre.append([frontfind0,samechre,i,bothscore])
			
	allsamechre=sorted(allsamechre,key=lambda x : x[-1], reverse=True)[:20]

	for allsamechre0 in allsamechre:

		frontfind0=allsamechre0[0]
		samechre=allsamechre0[1]		
		i=allsamechre0[2]

		chr0,std,s0,e0=frontfind0[0],frontfind0[1],frontfind0[2],frontfind0[3]

		if len(samechre)==0:

			continue
			
		
		else:
			
			find=max(1,find)
			
			match0=samechre[0]

			cutstart=max(0,min(s0,e0,match0[2],match0[3])-50000)

			cutend=max(s0,e0, match0[2],match0[3])+50000

			examscore=examposi(seqs,assemseqs,std,[chr0,e0,match0[2]],[contig0,codis[1],codis[2]],svfolder+name+'_exam'+str(i),lf-fn,le-en)


			if std=='+':

				cutmrna(seqs,chr0,std,cutstart,cutend,svfolder+name+'_retrorange'+str(i))

				rangestart=cutstart
		
				rangeend=cutend

			else:

				cutmrna(seqs,chr0,std,cutend,cutstart,svfolder+name+'_retrorange'+str(i))

				rangestart=cutend

				rangeend=cutstart

			os.system('makeblastdb -in {:s} -dbtype nucl -parse_seqids -out {:s}'.format(svfolder+name+'_retrorange'+str(i),svfolder+name+'_retrorange'+str(i)+'db'))

			score_dif=0
			
			realgap=[-1,-1,-1,-1,0.0,-1,-1,-1,100*rn/lr,-1]

			fullfind,full0,fullN=blastncheck('denova2/{:s}/'.format(inputfile)+name+'_fullseq',svfolder+name+'_retrorange'+str(i)+'db','denova2/{:s}/'.format(inputfile)+name+'_fullseqblastn'+str(i),50,50)

			fullfind=sorted([x for x in fullfind if x[-1]>30 and x[-3]<lf and x[-2]<le],key=lambda x: x[-1])

			if len(fullfind)==0:

				if std=='+':

					sd2=1

				else:

					sd2=-1


				score_dif,realgap=findinsertion(svfolder,inputfile,name,i,svfolder+name+'_retrorange'+str(i),lf,lfr,rn)

			
				if realgap[-2] != -1:

					realgap[0]=rangestart+sd2*realgap[0]
					realgap[1]=rangestart+sd2*realgap[1]


					realgap[2]=realgap[2]+fulls
					realgap[3]=realgap[3]+fulls

				else:

					realgap[-2]==100.0*fullN/full0


			else:


				if std=='+':

					sd1=1

				else:

					sd1=-1

				fullfind=fullfind[-1]			

				if fullfind[1]=='+':

					sd=1

				else:

					sd=-1


				if (std=='+' and fullfind[1]=='+') or (std=='-' and fullfind[1]=='-'):

					fullfind[1]='+'
					sd2=1					
					refstart=max(0,rangestart+sd1*fullfind[2]-lf)
					refend=rangestart+sd1*fullfind[3]+le


				else:

					fullfind[1]='-'
					sd2=-1
					refend=max(0,rangestart+sd1*fullfind[3]-le)
					refstart=rangestart+sd1*fullfind[2]+lf

				cutmrna(seqs,chr0,fullfind[1],refstart, refend,'denova2/{:s}/'.format(inputfile)+name+'_fullref'+str(i))

				score_dif,realgap=findinsertion(svfolder,inputfile,name,i,'denova2/{:s}/'.format(inputfile)+name+'_fullref'+str(i),lf,lfr,rn)

				if realgap[-2] != -1:
				
					realgap[0]=refstart+sd2*realgap[0]
					realgap[1]=refstart+sd2*realgap[1]


					realgap[2]=realgap[2]+fulls
					realgap[3]=realgap[3]+fulls

				else:
				
					realgap[-2]=100.0*fullN/full0
	
			realgap.insert(2, contig0)	
			
			retrofind=findretro(svfolder,svfolder+name+'_retrorange'+str(i),name,i,chr0,std,rangestart,rangeend,seqs)	

			out=[name,i,chr0]+retrofind[:5]+[frontfind0[6],match0[6]]+retrofind[5:]+[retrosize-retroN,lf-fn,le-en]+realgap+[score_dif,examscore]

			allscore.append(examscore)
			
			allout.append(out)

			if '_' not in chr0:
				mainscore.append(examscore)
			else:

				mainscore.append('a')

	lock1.acquire()


	if len(allout)==0:

		write=pd.DataFrame.from_records([[name,0,retrosize,lf0,le0,retroN,frontN,endN]]).to_csv('/home/q5476572/data/denova2/'+inputfile+'_unknown',mode='a',sep='\t',header=None,index=False)

		lock1.release()

		return -1


	write=pd.DataFrame.from_records(allout).to_csv('/home/q5476572/data/denova2/'+inputfile+'_allout',mode='a',sep='\t',header=None,index=False)

	best=max(allscore)

	if len([x for x in mainscore if type(x)!=type('a')])>0:

		bestmain=max([x for x in mainscore if type(x)!=type('a')])

		if bestmain >= 60:

			best=bestmain

	goodout=[x for x in allout if x[-1]>=max(60, 0.9*best)]


	if len(goodout)==0:

		write=pd.DataFrame.from_records([[name,0,retrosize,lf0,le0,retroN,frontN,endN]]).to_csv('/home/q5476572/data/denova2/'+inputfile+'_unknown',mode='a',sep='\t',header=None,index=False)

		lock1.release()

		return -1

	goodnega=[x for x in goodout if (x[6]*x[7]<=8000 or x[10]=='shift') and  (x[-8]*x[-3]<=8000)]

	goodposi=[x for x in goodout if (x[6]*x[7]>8000 and x[10]=='-') and (x[-3]*x[-8]>8000)]

	unknown=[x for x in goodout if x not in goodnega+goodposi]

	goodposimain=[x for x in goodposi if '_' not in x[2]]


	if len(goodposimain)>0:

		write=pd.DataFrame.from_records(goodposimain).to_csv('/home/q5476572/data/denova2/'+inputfile+'_posi',mode='a',sep='\t',header=None,index=False)

	else:
		if len(unknown)>0:

			write=pd.DataFrame.from_records(unknown).to_csv('/home/q5476572/data/denova2/'+inputfile+'_tocheck',mode='a',sep='\t',header=None,index=False)

		if len(goodposi)>0:

			write=pd.DataFrame.from_records(goodposi).to_csv('/home/q5476572/data/denova2/'+inputfile+'_posi',mode='a',sep='\t',header=None,index=False)

		if len(goodnega)>0:
			
			write=pd.DataFrame.from_records(goodnega).to_csv('/home/q5476572/data/denova2/'+inputfile+'_nega',mode='a',sep='\t',header=None,index=False)


	lock1.release()
		
	return 0



def runblastn2(arg):
	
	global seqs,assemseqs

	name,inputfile,contig,start,end=arg[0],arg[1],arg[2],arg[3],arg[4]
	
	svfolder='temp/{:s}/'.format(inputfile)
	
	rangefile='db/hg38'
	
	with open(svfolder+name+'_f1',mode='r') as f:

		read0=f.read()
		fronttitle=read0.split('\n')[0]
		frontseq=''.join(read0.split('\n')[1:])

	f.close()

	del read0

	with open(svfolder+name+'_e1',mode='r') as f:

		read0=f.read()
		endtitle=read0.split('\n')[0]
		endseq=''.join(read0.split('\n')[1:])

	f.close()

	del read0
	
	lf=len(frontseq)
	
	le=len(endseq)

	fn,en=frontseq.count('N'),endseq.count('N')
	
	with open('denova2/{:s}/'.format(inputfile)+name+'_fullseq',mode='r') as f:

		read0=f.read()
		fronttitle=read0.split('\n')[0]
		frontseq=''.join(read0.split('\n')[1:])

	f.close()
	
	with open(svfolder+name+'_retro',mode='r') as f:

		read0=f.read()
		retroseq=''.join(read0.split('\n')[1:])

	f.close()
	
	retrosize=len(retroseq)

	lfr=lf+retrosize
	
	retroN=retroseq.count('N')
	
	contig0='_'.join(fronttitle.split('_')[:-2])[1:]
	
	fulls=int(fronttitle.split('_')[-2])
	
	fulle=int(fronttitle.split('_')[-1])
	

	find,full0,fullN=blastncheck('denova2/{:s}/'.format(inputfile)+name+'_fullseq',rangefile,svfolder+name+'_fullout')


	find=sorted([x for x in find if x[-1]>30 and x[-3]<lf and x[-2]<le],key=lambda x : x[-1], reverse=True)[:20]

	allout=[]
	for i,find0 in enumerate(find):
		
		chr0,std,s0,e0=find0[0],find0[1],find0[2],find0[3]

		if std=='+':
			sd=1
			cutstart0=max(0,s0-50000)
			cutend0=e0+50000
		else:
			sd=-1
			cutstart0=s0+50000
			cutend0=max(0,e0-50000)

		cutstart=s0

		cutend=e0
	
		fullscore=examposi(seqs,assemseqs,std,[chr0,s0,e0],[contig0,fulls,fulle],'denova2/{:s}/'.format(inputfile)+name+'_'+str(i)+'_fullseqcheckout',lf,le)
		
		cutmrna(seqs,chr0,std,cutstart0,cutend0,'denova2/{:s}/'.format(inputfile)+name+'_fullref2'+str(i))

		score_dif,realgap=findinsertion(svfolder,inputfile,name,i,'denova2/{:s}/'.format(inputfile)+name+'_fullref2'+str(i),lf-find0[-2],lfr-find0[-2],retroN)
	
		os.system('makeblastdb -in {:s} -dbtype nucl -parse_seqids -out {:s}'.format('denova2/{:s}/'.format(inputfile)+name+'_fullref2'+str(i),'denova2/{:s}/'.format(inputfile)+name+'_fullref2'+str(i)+'db'))

		retrofind=findretro(svfolder,'denova2/{:s}/'.format(inputfile)+name+'_fullref2'+str(i),name,i,chr0,std,cutstart0,cutend0,seqs)

		if std=='+':
			sd=1
		else:
			sd=-1
	
		frontsize=realgap[0]
		backsize=abs(e0-s0)-realgap[1]


		if realgap[-2] != -1:
		
			realgap[0]=sd*realgap[0]+cutstart
			realgap[1]=sd*realgap[1]+cutstart
			realgap[2]=realgap[2]+fulls
			realgap[3]=realgap[3]+fulls
			
			
		else:
			
			realgap[-2]=100.0*fullN/full0
			
		realgap.insert(2, contig0)
			
		out=[name,i,chr0]+retrofind[:5]+[find0[-1],find0[-1]]+retrofind[5:]+[retrosize-retroN,frontsize, backsize]+realgap+[score_dif,fullscore]
		
		allout.append(out)
		
	
	if len(allout)==0:

		return -1

	lock1.acquire()

	write=pd.DataFrame.from_records(allout).to_csv('/home/q5476572/data/denova2/'+inputfile+'_allout',mode='a',sep='\t',header=None,index=False)

	lock1.release()



	goodanchor=len([x for x in [lf-fn,le-en] if x>=100])
	
	best=max([x[-1] for x in allout])

	main=[x for x in allout if '_' not in x[2]]

	if len(main)>0:

		bestmain=max([x[-1] for x in main])

		if bestmain>=30*goodanchor:

			bestmain=best

	goodout=[x for x in allout if x[-1]>=max(30*goodanchor, 0.9*best)]	

	if len(goodout)==0:

		return -1	

	nega=[x for x in goodout if (x[6]*x[7]<=8000 or x[10]=='shift') and  (x[-3]*x[-8]<=8000)]

	posi=[x for x in goodout if (x[6]*x[7]>8000 and x[10]=='-') and (x[-3]*x[-8]>8000)]

	unknown=[x for x in goodout if x not in nega+posi]

	posimain=[x for x in posi if '_' not in x[2]]

	

	lock1.acquire()

	if len(posimain)>0:

		write=pd.DataFrame.from_records(posimain).to_csv('/home/q5476572/data/denova2/'+inputfile+'_posi',mode='a',sep='\t',header=None,index=False)

	else:
		
		
		if len(unknown)>0:

			write=pd.DataFrame.from_records(unknown).to_csv('/home/q5476572/data/denova2/'+inputfile+'_tocheck',mode='a',sep='\t',header=None,index=False)

		if len(posi)>0:

			write=pd.DataFrame.from_records(posi).to_csv('/home/q5476572/data/denova2/'+inputfile+'_posi',mode='a',sep='\t',header=None,index=False)

		if len(nega)>0:
			
			write=pd.DataFrame.from_records(nega).to_csv('/home/q5476572/data/denova2/'+inputfile+'_nega',mode='a',sep='\t',header=None,index=False)
	

	lock1.release()
		
	return 0


	


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

	t1=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_notfind' ,sep=',',header=None)

	sort1=list(t1[0])

except:

	sort1=[]

try:
	t2=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_find' ,sep=',',header=None)

	sort2=list(t2[0])

except:
		
	sort2=[]

try:
	t3=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_singlenega' ,sep=',',header=None)

	sort3=list(set(list(t3[0])))

except:
	sort3=[]


try:
	t4=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_singleposi' ,sep=',',header=None)

	sort4=list(set(list(t4[0])))

except:

	sort4=[]

t0=pd.read_csv(inputfile,sep=',')

presort=list(t0['sort'])

presort=[i for i in xrange(len(t0)) if presort[i] in sort1+sort2+sort3 and presort[i] not in sort4]

t0=t0.iloc[presort]

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

#os.system('python cutmrnas.py -f {:s} -s {:d} -c'.format(inputfile, 30000))

for file0 in [x for x in  os.listdir('chroms') if '.fa' in x]:

	with open('chroms/'+file0, mode='r') as f:
		read=f.read()
	f.close()
	
	seqs[read.split('\n')[0][1:]]=''.join(read.split('\n')[1:])

	del read

if assem in os.listdir(assemfolder):
	
	with open(assemfolder+assem, mode='r') as f:
		read=f.read().split('>')[1:]
	f.close()
	
	for read0 in read:
		
		assemseqs[read0.split('\n')[0].split(' ')[0]]=''.join(read0.split('\n')[1:])

	del read


args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]), 30000] for i in xrange(len(t0))]
test=0
if test==0:
	m0=4
	p=mul.Pool(processes=m0)
	left=p.map_async(runblastn,args)
	p.close()
	p.join()


	left=left.get()	
elif test==1:
	print 'a'
	map(runblastn,args)


try:	
	t1=pd.read_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_tocheck' ,sep='\t',header=None)
	sort1=list(t1[0])
except:
	sort1=[]

	
try:	
	t2=pd.read_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_unknown' ,sep='\t',header=None)
	sort2=list(t2[0])
except:
	sort2=[]

unknown=list(set([str(x) for x in sort2]))

args = [[str(sort[i]),inputfile.split('/')[-1],str(contig[i]),int(assemble_start[i]),int(assemble_end[i]), 30000] for i in xrange(len(t0)) if str(sort[i]) in unknown]

if test==2:
	left=map(runblastn2,args)

m0=4
p=mul.Pool(processes=m0)
left=p.map_async(runblastn2,args)
p.close()
p.join()

left=left.get()	
leftsort=[args[i][0] for i in xrange(len(left)) if left[i]<0]

index1=[i for i in xrange(len(sort1)) if sort1[i] in leftsort]
index2=[i for i in xrange(len(sort2)) if sort2[i] in leftsort]

write=t2.iloc[index2].to_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_unknown',mode='w',sep='\t',header=None,index=False)




t1=pd.read_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_nega' ,sep='\t',header=None)

sort1=list(t1[0])

t2=pd.read_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_posi' ,sep='\t',header=None)

sort2=list(t2[0])
chrom2=list(t2[2])
samechrom=cl.defaultdict(list)

for i,x in enumerate(sort2):

	samechrom[x].append(chrom2[i])


alter=[]
for x in samechrom.keys():

	if len([a for a in samechrom[x] if '_' not in a])==0:

		alter.append(x) 

isalterposi=[i for i  in xrange(len(sort2)) if sort2[i] in alter]

try:
	t3=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_alter' ,sep=',',header=None)

	sort3=list(t3[0])
except:

	sort3=[]

try:
	t4=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_singlealter' ,sep=',',header=None)

	sort4=list(set(list(t4[0])))

except:

	sort4=[]

try:
        t5=pd.read_csv('/home/q5476572/data/denova/'+inputfile.split('/')[-1]+'_singleposi' ,sep=',',header=None)

        sort5=list(set(list(t5[0])))

except:

        sort5=[]



notalter=[i for i in xrange(len(sort1)) if sort1[i] not in sort2+sort3+sort4]


isalter=[i for i in xrange(len(sort3)) if (sort3[i] not in sort2 or sort3[i] in alter) and sort3[i] not in sort5]

isalter2=[i for i in xrange(len(sort4)) if (sort4[i] not in sort2 or sort4[i] in alter) and sort4[i] not in sort5]












if len(t1) >0:

	write=t1.iloc[notalter].to_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_denovaindex',mode='w',sep='\t',header=None,index=False)

if len(t2)>0:

	write=t2.iloc[isalterposi].to_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_alterindex',mode='w',sep='\t',header=None,index=False)

if len(t3)>0:

	write=t3.iloc[isalter].to_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_alterindex2',mode='w',sep='\t',header=None,index=False)	

	write=t4.iloc[isalter2].to_csv('/home/q5476572/data/denova2/'+inputfile.split('/')[-1]+'_alterindex2',mode='a',sep='\t',header=None,index=False) 


