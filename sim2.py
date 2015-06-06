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


from sys import argv
from scarab_classes import *
from graphics import *
from time import sleep
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D


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
    hp=0
    for k in habitat_probs:
        if habitat_probs[k]==h_max:
            hp=k     
    habitat_area=[]    
    for a in areas.values():
        if a['habitat']==hp:
            habitat_area+=a['area']
       
    
    p=world.random_positions(n=n,area=habitat_area)     
    world.create_agents(agent_type=Beetle, n=n, pos=p, **args)



grid_size=(200,200)
#traps_file='/home/fabio/Desktop/monografia/spillover/sampling-2'
#run=1
#days=2


#sampling_design=traps_file.split('/')
#sampling_design=sampling_design[len(sampling_design)-1]


#habitat_layer=load_habitats('/home/fabio/Desktop/monografia/spillover/land-2')
surface=np.zeros((2,grid_size[0],grid_size[1]))
#fragment1=[k for k in habitat_layer.keys() if habitat_layer[k]==54]
#fragment2=[k for k in habitat_layer.keys() if habitat_layer[k]==84]
#matrix=[k for k in habitat_layer.keys() if habitat_layer[k]==0]

areas={}
areas[0]={'habitat':1, 'area':[(x,y) for x in xrange(grid_size[0]) for y in xrange(grid_size[1])]}
#areas[1]={'habitat':1,'area':fragment1}
#areas[2]={'habitat':1,'area':fragment2}
#areas[3]={'habitat':0,'area':matrix}


habitat_array=fill_habitat_array(areas, size=(200,200))

#surface[0]=habitat_array
grid=Rectangular_Grid(x_max=200,y_max=200,areas=areas,ndim=2,dim_names=("Habitat","Cues"))
grid.load_surface(surface,('Habitat','Cues'))




world=World(grid,steps_in_day=1)


#trap_args={'grid':grid,'world':world,
#'dim':'Cues','reach': 8,'radius':5.0}

#trap_pos=load_trap_pos(traps_file)
#trap_pos=[(100,100)]
#world.create_agents(agent_type=Trap, n=len(trap_pos), pos=trap_pos, **trap_args)

#dung_pos=world.random_positions(n=5,area=areas[1]['area'])
#dung_pos+=world.random_positions(n=5,area=areas[2]['area'])
#dung_pos+=world.random_positions(n=15,area=areas[3]['area'])
dung_pos=world.random_positions(n=15,area=areas[0]['area'])

#dung_pos=[(100,100),(50,50),(100,90),(38,50)]

dung_args={'grid':grid,'world':world,
'dim':'Cues','reach': 10,'radius':5.0,'amount':2,'max_age':5}

world.create_agents(agent_type=Dung, n=len(dung_pos), pos=dung_pos, **dung_args)

#load_spp_par('/home/fabio/Desktop/monografia/spillover/spp_par',world)
active_days=range(0,180)+range(366,366+180)
#active_days=range(180,730)
breeding_days=range(90,180)+range(366+90,366+180)
sp1_args={'grid':grid,'sp':'sp.1','sex':'F',
                'world':world,'habitat_probs':{0:0.999,1:0.001},
                'dist_par':{'mean':5,'sd':1,'max_dist':3000},
                'energy_par':{'hungry':15, 'max_energy':30,'initial_energy':25,
                'rest':0.02, 'move':0.5, 'breed':4}, 'min_age':60,'age':0, 'max_age':1000,
                'activity':{0:False,1:True,2:True,3:True,4:True,5:True,6:True,
                7:True,8:True,9:True,10:True,11:True,12:True,13:True,14:True,
                15:True,16:True,17:True,18:True,19:True,20:True,21:False,22:False,
                23:False,24:False},'active_days':active_days,'breeding_days':breeding_days, 'perception_radius':5}

p=world.random_positions(n=100,area=areas[0]['area'])
#p=[(93,93),(100,95),(40,50),(70,50)]
world.create_agents(agent_type=Beetle, n=len(p), pos=p, **sp1_args)
for b in Beetle.Instances.values():
    b.sex=np.random.choice(["F","M"])
    b.age=np.random.randint(365)
    b.energy=np.random.randint(10,30)
#b1=Beetle.Instances.values()[0]
#p=b1.grid.circle(center=b1.what_cell(b1.position),radius=b1.perception_radius)
#partner=b1.suitable_partner(p)
#b1.breed(partner)


N=[]
F1=[]
F2=[]
#####GUI#####
def main():
    
        
    win = GraphWin("Sim 2",grid.x_max, grid.y_max)
    for x in xrange(grid.x_max):
        for y in xrange(grid.y_max):
            if grid.surface[0][x,y]==1:
                win.plot(x,y,"green")
            
    
            
        
    
    for t in Trap.Instances.values():    
        c = Circle(Point(t.position[0],t.position[1]), 1)
        c.setFill('black')
        c.draw(win)
    '''     
    Dung_dots={}    
    for d in Dung.Instances.values():    
        c = Circle(Point(d.position[0],d.position[1]), 1)
        c.setFill('red')
        c.setOutline('red')
        c.draw(win)
        Dung_dots[d.id]=c
    ''' 
    Beetle_dots={}       
    for b in Beetle.Instances.values():    
        c = Circle(Point(int(b.position[0]),int(b.position[1])), 1)
        c.setFill('blue')
        c.setOutline('blue')
        c.draw(win)
        #print b.position
        #c.move_to(b.position[0],b.position[1])
        Beetle_dots[b.id]=c
    
      
    #world.step=25   
    Dung_dots={}    
    for i in xrange(365*2):
        print world.day, Beetle.PopulationSize(), world.report_functions()[0].values()[0]
        #dung_pos=world.random_positions(n=1,area=areas[1]['area'])
        #dung_pos+=world.random_positions(n=1,area=areas[2]['area'])
        #dung_pos+=world.random_positions(n=5,area=areas[3]['area'])
        dung_pos=world.random_positions(n=3,area=areas[0]['area'])
        dung_args={'grid':grid,'world':world,
        'dim':'Cues','reach': 10,'radius':5.0,'amount':1.8,'max_age':7}

        world.create_agents(agent_type=Dung, n=len(dung_pos), pos=dung_pos, **dung_args)
        
        for d in Dung.Instances.values():    
            if d.id not in Dung_dots.keys(): 
                c = Circle(Point(d.position[0],d.position[1]), 1)
                c.setFill('red')
                c.setOutline('red')
                c.draw(win)
                Dung_dots[d.id]=c
    
        for b in Beetle.Instances.values():    
            if b.id not in Beetle_dots.keys():
                c = Circle(Point(d.position[0],d.position[1]), 1)
                c.setFill('blue')
                c.setOutline('blue')
                c.draw(win)
                Beetle_dots[b.id]=c 
        
        shuffled_instances=Beetle.Instances.values()
        np.random.shuffle(shuffled_instances)        
        for b in shuffled_instances: 
            
            #b.orientation=angle_cells(b.position,(100,100))%360
            #b.move(1, angle=b.orientation)
            #sleep(.05)
            b.action()
            #print b.id, b.position, b.energy
            #sleep(0.5)
            #print b.energy
            #b.wiggle(1,sd=80)
            if b.age<b.min_age:
                Beetle_dots[b.id].setFill('yellow')
                Beetle_dots[b.id].setOutline('yellow')
            else:
                Beetle_dots[b.id].setFill('blue')
                Beetle_dots[b.id].setOutline('blue')    
            
            Beetle_dots[b.id].move_to(int(b.position[0]),int(b.position[1]))
        
        for bd in Beetle_dots.keys():
            if bd not in Beetle.Instances.keys(): Beetle_dots[bd].undraw()    
        for dd in Dung_dots.keys():
            if dd not in Dung.Instances.keys(): Dung_dots[dd].undraw()    

        
        n=Beetle.PopulationSize()
        N.append(n)
        F1.append(world.report_functions()[0].values()[0])
        F2.append(world.report_functions()[0].values()[1])
        nr=ceil(n*0.01)
        if nr<n and np.random.random()<0.05: Beetle.RemoveRandomly(nr)         
        world.check_traps()
        world.increment_time()
        
        
    

    '''
    
    world.step=25
    b1.orientation=0
    for i in xrange(5):
        b1.wiggle(1)
        print b1.orientation
        for i in xrange(100):
            b1.move(1,b1.orientation)
            Beetle_dots[b1.id].move_to(int(b.position[0]),int(b.position[1]))
            #world.increment_time()
            sleep(.2)
         
   
    
    b1=Beetle.Instances.values()[0]
    
    
    world.step=25
    for i in xrange(5000):
            b1.wiggle(1)
            Beetle_dots[b1.id].move_to(int(b1.position[0]),int(b1.position[1]))
            #world.increment_time()
            
         
   ''' 
    
    







#####################
