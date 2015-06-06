'''
This file is part of Scarab_IBM.
    
    Scarab-IBM provides a plataform for developing individual-based models
    for Dung Beetles (Coleoptera:Scarabaeidae:Scarabaeinae)
    
    Copyright (C) 2015  Fabio Soares Frazao
    
    Scarab_IBM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Scarab-IBM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Scarab-IBM.  If not, see <http://www.gnu.org/licenses/>.
'''



import numpy as np
from random import random
from math import ceil,floor,atan2,degrees,sin,cos,radians,tan,pi
import csv

#CONVENIENCE FUNCTIONS

#Cauchy generator 
def cauchy(location, scale):
 
    # Start with a uniform random sample from the open interval (0, 1).
    # But random() returns a sample from the half-open interval [0, 1).
    # In the unlikely event that random() returns 0, try again.
     
    p = 0.0
    while p == 0.0:
        p = random()
         
    return location + scale*tan(pi*(p - 0.5))

#FILE PROCESSING

def array2csv(array,filename):
    np.savetxt(filename,array,delimiter=',')
    
        
def csv2array(filename):
    return np.genfromtxt(filename,delimiter=',')    
    
def WriteReport(filename,dic,sep,rownames=['dung','pedoturbation','soil_revolved']):
    with open(filename,'a') as f:
        f.write(rownames[0]+sep+sep.join([n for n in rownames[1:]])+'\n')
        for k in dic.keys():
            f.write(str(k) + sep +sep.join([str(x) for x in dic[k].values()]) + '\n')
    

def load_habitats(filename):
    with open(filename,'r') as f:
        cells=f.readline().split()
    habitats={ (abs(int(cells[c])),abs(int(cells[c+1]))):int(cells[c+2]) for c in xrange(0,len(cells),3)}    
    return habitats 

    
        
    
def fill_habitat_array(areas,size):
    array=np.zeros(size)
    for a in areas.values():
        for cell in a['area']: array[cell]=int(a['habitat'])
     
    return array
    
    
def load_trap_pos(filename):
    pos=[]
    with open(filename,'r') as f:
           for line in f:
            line_content=line.split()
            pos.append((abs(int(line_content[0])),abs(int(line_content[1]))))     
    return pos
    
def ask(agents, methodname, *args, **kwargs):
	f = methodcaller(methodname, *args, **kwargs)
	for agent in agents:
			f(agent)
			
def euclid_dist(x1,x2,y1,y2):
    return np.sqrt( ((x1-x2)**2) + ((y1-y2)**2) )			


def angle_cells((x0,y0),(x1,y1)):
    deltax=x0-x1
    deltay=y0-y1
    
    return degrees(atan2(deltay,deltax))
    
#Movement			
def quad(angle):
    if angle>=360 or angle<=0:
        angle=360
    
    quads={1:(0,90),2:(90,180),3:(180,270),4:(270,360)}
    for q in quads.keys():
        min_a,max_a=quads[q]
        if angle>min_a and angle<=max_a: return q
    
def new_pos(x0,y0,d,angle):
       
    return ((x0+cos(radians(angle))*d,y0+sin(radians(angle))*d))
    

def load_spp_par(filename,world):
    #sp_par=[]
    with open(filename,'r') as f:
        for line in f:
            sp_par=eval(line)
            distribute_beetles(sp_par,world)
    
def distribute_beetles(sp_par,world):
    areas=world.grid.areas
    n=sp_par['n']
    args=sp_par['args']
    args['sex']=np.random.choice(["M","F"])
    habitat_probs=args['habitat_probs']
    h_max=max(habitat_probs.values())
    hp=[]
    for k in habitat_probs:
        if habitat_probs[k]==h_max:
            hp.append(k)     
    habitat_area=[]    
    for a in areas.values():
        if a['habitat'] in hp:
            habitat_area+=a['area']
       
    
    p=world.random_positions(n=n,area=habitat_area)     
    world.create_agents(agent_type=Beetle, n=n, pos=p, **args)



