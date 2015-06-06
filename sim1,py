'''

    This program uses Scarab_IBM to implement Simulation 1.
    
    Copyright (C) 2015  Fabio Soares Frazao

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#import numpypy
from sys import argv
from scarab_classes import *
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D




#'''
traps_file='/home/fabio/Desktop/monografia/spillover/sampling-2'
run=1
days=2
#'''

'''
traps_file=argv[1]
run=argv[2]
days=argv[3]
#'''
print str(days) + str(run)

sampling_design=traps_file.split('/')
sampling_design=sampling_design[len(sampling_design)-1]


habitat_layer=load_habitats('/home/fabio/Desktop/monografia/spillover/land-2')
surface=np.zeros((2,200,200))
fragment1=[k for k in habitat_layer.keys() if habitat_layer[k]==54]
fragment2=[k for k in habitat_layer.keys() if habitat_layer[k]==84]
matrix=[k for k in habitat_layer.keys() if habitat_layer[k]==0]

areas={}
areas[1]={'habitat':1,'area':fragment1}
areas[2]={'habitat':1,'area':fragment2}
areas[3]={'habitat':0,'area':matrix}


habitat_array=fill_habitat_array(areas, size=(200,200))

surface[0]=habitat_array
grid=Rectangular_Grid(x_max=200,y_max=200,areas=areas,ndim=2,dim_names=("Habitat","Cues"))
grid.load_surface(surface,('Habitat','Cues'))




world=World(grid,steps_in_day=500)


trap_args={'grid':grid,'world':world,
'dim':'Cues','reach': 10,'radius':0.7}

trap_pos=load_trap_pos(traps_file)
world.create_agents(agent_type=Trap, n=len(trap_pos), pos=trap_pos, **trap_args)


dung_pos=world.random_positions(n=10,area=areas[1]['area'])
dung_pos+=world.random_positions(n=10,area=areas[2]['area'])
dung_pos+=world.random_positions(n=40,area=areas[3]['area'])

dung_args={'grid':grid,'world':world,
'dim':'Cues','reach': 10,'radius':0.3,'amount':20,'max_age':6}

world.create_agents(agent_type=Dung, n=len(dung_pos), pos=dung_pos, **dung_args)


np.random.seed(10)
load_spp_par('/home/fabio/Desktop/monografia/spillover/spp_par',world)
#np.random.seed()


community=world.report_community()
rownames=Beetle.CountBeetles(Beetle.Instances.values()).keys()
rownames.insert(0,'Area')
WriteReport('/home/fabio/Desktop/monografia/spillover/results/community_' + sampling_design +'_'+ str(days)+ '_' +  str(run) + '.txt',
dic=community,sep=';',rownames=rownames)


for i in xrange(200):
    for b in Beetle.Instances.values(): b.action2()
    world.increment_time()



for i in xrange(world.steps_in_day*int(days)):
    for b in Beetle.Instances.values(): b.action2()
    world.increment_time()
    world.check_traps()
    
samples=world.report_samples()
rownames=Beetle.CountBeetles(Beetle.Instances.values()).keys()
rownames.insert(0,'Trap')
WriteReport('/home/fabio/Desktop/monografia/spillover/results/sample_' + sampling_design + '_' + str(days) + '_' +  str(run) + '.txt',
dic=samples,sep=';',rownames=rownames)
    
