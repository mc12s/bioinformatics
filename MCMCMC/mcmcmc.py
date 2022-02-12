#!/usr/bin/env/python
from random import uniform
from random import randint
from math import exp
import sys

step=int(sys.argv[1])
N=int(sys.argv[2])
m=int(sys.argv[3])

def swap(a,b):
	R=(feval[b][1]**Tinv[a]*feval[a][1]**Tinv[b])/(feval[b][1]**Tinv[b]*feval[a][1]**Tinv[a])
	#print "R=",R	
	if uniform(0.0,1.0)<R:
	#	print "Swapped ",a," and ",b	
		return True
	#else:
	#	print "No Swap"

def f(x,y):
	return exp((1/10.)*(-x**2-y**2)) + exp((1/200.)*(-(x-40.)**2-(y-40.)**2)) + .5*exp((1/20.)*(-(x-20.)**2-(y-20.)**2))

def accept(feval_try,feval_old,chain):
	if (uniform(0.0,1.0)<(feval_try/feval_old)**Tinv[chain]):
		return True
	
#Initialize heated chains
Tinv=[0.0 for i in range(m)]
Tdiv=1.0/m + .001
for i in xrange(0,m):
	Tinv[i]=1.-Tdiv*i

#Initialize positions for m chains at 0,0
pos_try=[[0.0 for j in range(2)] for i in range(m)]
pos_old=[[0.0 for j in range(2)] for i in range(m)]
feval=[[0.0 for j in range(2)] for i in range(m)]

index=0
#print Tinv
for i in xrange(1,N):
	accepted=0
	for chain in xrange(0,m):				#test acceptance for all chains
		r1,r2=uniform(-1,1),uniform(-1,1)
		pos_try[chain]=[pos_old[chain][0]+(r1*step),pos_old[chain][1]+(r2*step)]
		feval[chain]=[f(pos_try[chain][0],pos_try[chain][1]),f(pos_old[chain][0],pos_old[chain][1])]
		if accept(feval[chain][0],feval[chain][1],chain):
			if chain==0:
				accepted=1
			pos_old[chain]=pos_try[chain]
			feval[chain][1]=feval[chain][0]

	if (i%10000==0):
		#print i
		a,b=randint(0,m-1),randint(0,m-1)
		while a==b & m>1:
			b=randint(0,m-1)
		if swap(a,b):									#if swap then switch heats and index
			#print "Tinv(before)",Tinv		
			Tinv_temp=Tinv[a]
			Tinv[a]=Tinv[b]
			if Tinv[a]==1.0:
				index=a;
			Tinv[b]=Tinv_temp	
			if Tinv[b]==1.0:
				index=b;
			#print "Tinv(after)",Tinv
	print pos_old[index][0],pos_old[index][1],accepted;
