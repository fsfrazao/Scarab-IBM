

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

import scarab_helper

class Rectangular_Grid(object):

    def __init__(self,x_max,y_max,areas,ndim=1,dim_names=None):
        self.surface=np.zeros((ndim,x_max, y_max))
        self.x_max=x_max
        self.y_max=y_max
        self.dim_names={}
        self.areas=areas
        for i in xrange(0,ndim): self.dim_names[dim_names[i]]=i   

    
    def load_surface(self,surface,dim_names):
        self.surface=surface
        if len(surface.shape)==3:
            ndim=surface.shape[0]
            self.x_max=surface.shape[1]
            self.y_max=surface.shape[2]
            self.dim_names={}
            for i in xrange(0,ndim): self.dim_names[dim_names[i]]=i
    
        else:
            self.x_max=surface.shape[0]
            self.y_max=surface.shape[1]
            
        
    def value(self,dim,cell):
        return self.surface[self.dim_names[dim]][cell]
        
    def change_value(self,dim,cell,new_value):
        self.surface[self.dim_names[dim]][cell]=new_value
    
    
    def hood(self,cell,radius=1,remove_center=True):    
        hood=[(x,y) for x in xrange(cell[0]-radius,cell[0]+radius+1) for y in xrange(cell[1]-radius,cell[1]+radius+1) if self.in_bounds((x,y))]
        if remove_center: hood.remove(cell)
        return hood
        
    def circle(self,center,radius):

        xmax=center[0]+radius+1
        if xmax>self.x_max: xmax=self.x_max
        ymax=center[1]+radius+1
        if ymax>self.y_max: ymax=self.y_max
        xmin=center[0]-radius-1
        if xmin<0: xmin=0        
        ymin=center[1]-radius-1
        if ymin<0: ymin=0



        pixels_in_circle=[(x,y) for x in range(xmin,xmax) for y in range(ymin,ymax) if euclid_dist(x1=x,x2=center[0],y1=y,y2=center[1])<radius ]

        return pixels_in_circle    

    def set_plume(self,center,reach):

        for i in xrange(reach,0,-1):
            c=self.circle(center,i)
            for p in c:self.surface[self.dim_names["Cues"]][p]+=reach+1-i
            
    def remove_plume(self,center,reach):
        for i in xrange(reach,0,-1):
            c=self.circle(center,i)
            for p in c:self.surface[self.dim_names["Cues"]][p]-=reach+1-i
            
    def in_bounds(self,(x,y)):
            return x>=0 and x<self.x_max and y>=0 and y<self.y_max
            
    def adjust_to_bounds(self,(x,y)):
        if new_pos[0]<0: x=0
        if new_pos[1]<0: y=0
        if new_pos[0]>self.grid.x_max: x=self.grid.x_max
        if new_pos[1]>self.grid.y_max: y=self.grid.y_max
        
        return (x,y)
            
            
    def cell_with_value(self,cells,value,dim):
        return [c for c in cells if self.surface[self.dim_names[dim]][c]==value]
           


class World(object):
    
    """
    Beetles=[]
    Traps=[]
    Dungs=[]
    """
    
    
    
    def __init__(self,grid,steps_in_day=1440):
        self.grid=grid
        self.day=0
        self.step=0
        self.steps_in_day=steps_in_day
        self.Beetles={}
        self.Traps={}
        self.Dungs={}
        
    
    def random_positions(self,n,area):
        pos=[]
        for i in xrange(n):
            pos.append(area[np.random.randint(len(area))])
        return pos

    def create_agents(self,agent_type,n,pos, **kwargs):
        for i in xrange(n):
            agent_type(position=pos[i],**kwargs)
    
            
    def increment_time(self):
        self.step+=1
        if self.step==self.steps_in_day:
            self.step=0
            self.day+=1
            for b in Beetle.Instances.values(): b.distance_covered=0
            Dung.IncreaseAge()
            Dung.RemoveOld()
            Beetle.IncreaseAge()
            Beetle.RemoveStarving()           
            Beetle.RemoveOld()
            
            
            
            
    def report_functions(self):
        func={}
        func_names=Dung.Functions.values()[0].keys()
        
        for area in self.grid.areas:
            func[area]={f:0 for f in func_names}
            functions=Dung.FunctionsByDungSet(Dung.DungsInArea(self.grid.areas[area]['area']))
            for f,n in functions.items(): func[area][f]=n
            
        return func
    
    def spp_list(self):
        return Beetle.CountBeetles(Beetle.Instances.values()).keys()
    
    def report_community(self):
        community={}
        spp=self.spp_list()
        
        for area in self.grid.areas:
            community[area]={sp:0 for sp in spp}
            beetles=Beetle.CountBeetles(Beetle.BeetlesInArea(self.grid.areas[area]['area']))
            for sp,n in beetles.items(): community[area][sp]=n
            
        return community
        
    def report_samples(self):
        sample={}
        spp=self.spp_list()
        
        for t in Trap.Instances:
            sample[Trap.Instances[t].id]={sp:0 for sp in spp}
            beetles=Trap.Instances[t].captured_beetles
            for sp,n in beetles.items(): sample[t][sp]=n
            
        return sample
        
    
        
      
        
    def check_traps(self):
        for t in Trap.Instances.values():
            beetles=Beetle.BeetlesInArea([t.position])
            for b in beetles:
                if euclid_dist(x1=t.position[0],x2=b.position[0],y1=t.position[1],y2=b.position[1])<t.radius:
                    t.add_beetle(b)

        


class Agent(object):
        
    def __init__(self,position,grid,world):
        self.position=position
        self.grid=grid
        self.world=world
        #self.grid.space[self.position]=1
        
class Trap(Agent):
        
    ID=0
    Instances={}
    
    @classmethod
    def IncrementID(self):
        self.ID+=1
        
        
    
        
    def __init__(self,position,grid,world,dim,reach,radius):
        super(Trap,self).__init__(position,grid,world)
        self.reach=reach
        self.dim="Cues"
        self.radius=radius
        self.captured_beetles={}
        self.grid.set_plume(center=position,reach=reach)
        self.world=world
        self.id=Trap.ID
        self.IncrementID()
        self.Instances[self.id]=self
        
        
    def add_beetle(self,beetle):
        if not self.captured_beetles.has_key(beetle.sp):
            self.captured_beetles[beetle.sp]=1
        else:
            self.captured_beetles[beetle.sp]+=1
        
        Beetle.Instances.pop(beetle.id)


class Dung(Trap):
        
    Functions={}
    ID=0
    Instances={}
    spp={}
    
    
    
    @classmethod
    def IncrementID(self):
        self.ID+=1
        
        
    @classmethod
    def IncreaseAge(self):
        for b in self.Instances.values():
            b.age+=1
    
    @classmethod
    def RemoveOld(self):
        for d in self.Instances.values():
            if d.age>=d.max_age:
                d.grid.remove_plume(center=d.position,reach=d.reach)
                self.Instances.pop(d.id)
                
    @classmethod
    def DungsInArea(self,area):
        return [d for d in self.Instances.values() if d.position in area]
        
    @classmethod
    def FunctionsByDungSet(self,dungs):
        ks=Dung.Functions.values()[0].keys()
        func={k:0 for k in ks}
        for k in ks:
            for d in Dung.Functions.values():
                func[k]+=d[k]
            
        return func    
                    

      
    def __init__(self,position,grid,world,dim,reach,radius,amount,max_age):
        super(Trap,self).__init__(position,grid,world)
        self.reach=reach
        self.dim="Cues"
        self.radius=radius
        self.grid.set_plume(center=position,reach=reach)
        self.amount=amount
        self.pedoturbation=0
        self.incorporated_dung=0
        self.world=world
        self.id=Dung.ID
        self.age=0
        self.max_age=max_age
        self.Instances[self.id]=self
        Dung.Functions[Dung.ID]={"pedoturbation":0,"incorporated_dung":0}
        self.IncrementID()
        
            
            
    
        
class Beetle(Agent):

    ID=0    
    Instances={}
    Community={}
    
    @classmethod
    def BeetlesInArea(self,area):
        return [b for b in Beetle.Instances.values() if b.what_cell(b.position) in area]
    
    @classmethod
    def CountBeetles(self,beetles):
        spp={}
        for b in beetles:
            if not spp.has_key(b.sp):
                spp[b.sp]=1
            else:
                spp[b.sp]+=1
        return spp        
    
    
    
    
    @classmethod
    def IncrementID(self):
        self.ID+=1
        
    @classmethod
    def PopulationSize(self):
        return len(self.Instances)
        
    @classmethod
    def IncreaseAge(self):
        for b in self.Instances.values():
            b.age+=1
            
    @classmethod
    def RemoveStarving(self):
        for b in self.Instances.values():
            if b.energy<=0:
                self.Instances.pop(b.id)

    
    @classmethod
    def RemoveOld(self):
        for b in self.Instances.values():
            if b.age>=b.max_age:
                self.Instances.pop(b.id)

    @classmethod
    def RemoveRandomly(self,n):
        ind=np.random.randint(low=0,high=len(self.Instances), size=n)
        for i in ind: self.Instances.pop(self.Instances.keys()[i])

    def __init__(self,position,grid,world,sp,sex,activity, active_days,breeding_days,habitat_probs,energy_par,dist_par,max_age,min_age,age, perception_radius):
        """dist_par: dictionary with mean and stardard deviation 
                     of the distribution from which distance will be
                      withdrawn. Ex {mean:0.7,sd:0.08, max_dist}
           energy_par: dictionary with energy parameters
                       {hungry:15, max_energy:30, initial_energy:20,
                        rest:0.1, move:0.2, breed:5}            
        
        """
        super(Beetle,self).__init__(position,grid,world)
        self.sp=sp
        self.sex=sex
        self.active_days=active_days
        self.breeding_days=breeding_days
        self.age=age
        self.perception_radius=perception_radius
        self.max_age=max_age
        self.min_age=min_age
        self.activity=activity
        self.max_dist=dist_par['max_dist']
        self.habitat_probs=habitat_probs
        self.distance_covered=0
        self.dist_par=dist_par
        self.world=world
        self.id=Beetle.ID
        self.IncrementID()
        self.Instances[self.id]=self    
        self.track=[]
        self.rest_time=self.resting_min(self.activity)
        self.orientation=np.random.choice(range(1,360),1)
        self.energy_par=energy_par
        self.energy=energy_par['initial_energy']
        
        
    def resting_min(self,dic):
        steps_in_hour=self.world.steps_in_day/24
        
        resting_min={k2:range(v2[0],v2[1])  for k2,v2 in {k1:(k1*steps_in_hour,(k1+1)*steps_in_hour)
        for k1 in dic.keys() if dic[k1]==False}.items()}    
        resting_minutes=[]
        for min_per_hour in resting_min.values():
            for mins in min_per_hour:
                #if mins<self.world.steps_in_day:
                resting_minutes.append(mins)
        resting_minutes.sort()    
        
        return resting_minutes
        
        
       
    def what_cell(self,(x,y)):
        return (int(floor(x)),int(floor(y)))
    
        
    def update_position(self,new_pos):
        x,y=new_pos
        if not self.grid.in_bounds(new_pos):
            if new_pos[0]<0: x=0+(0-new_pos[0])
            if new_pos[1]<0: y=0+(0-new_pos[1])
            if new_pos[0]>self.grid.x_max: x=self.grid.x_max-(new_pos[0]-self.grid.x_max)
            if new_pos[1]>self.grid.y_max: y=self.grid.y_max-(new_pos[1]-self.grid.y_max)
            self.orientation+=180
            self.orientation=self.orientation%360
            #self.wiggle(0)
        
        self.position=(x,y)
    
    def breed(self, partner,cells):
        sp_args={'grid':self.grid,'sp':self.sp,'sex':np.random.choice(['F','M']),
                'world':self.world,'habitat_probs':self.habitat_probs,
                'dist_par':self.dist_par,'energy_par':self.energy_par,'min_age':self.min_age,'age':0, 'max_age':self.max_age,
                'activity':self.activity,'active_days':self.active_days,'breeding_days':self.breeding_days,'perception_radius':self.perception_radius}

        
        if len(cells)>1:
            p=self.world.random_positions(n=1,area=cells)
            self.world.create_agents(agent_type=Beetle, n=len(p), pos=p, **sp_args)
            self.energy-=self.energy_par['breed']
            partner.energy-=partner.energy_par['breed']

        
    def move(self, distance, angle):
        self.orientation=angle%360
        npos=new_pos(self.position[0],self.position[1],distance,self.orientation)
        self.update_position(npos)
        self.distance_covered+=distance
        #self.track.append(self.position)
        
    def move_towards(self,d,aim_position):
        self.move(d,(angle_cells(aim_position,self.position))%360)
    
    def cell_with_higher_value(self, cells):
        if len(cells) > 1:
            hv_cells = {cell: value for cell, value in cells.items() if value == max(cells.values())}
            return hv_cells.items()[np.random.choice(range(len(hv_cells)))]
        else:
            return self.what_cell(self.position), 0.0
        
    def which_cell(self, cells,d):
        
        if len(cells)>1:        
            #hv_cells={cell:value for cell,value in cells.items() if value==max(cells.values())}
            
            x,y=self.cell_center(self.position)
            
            npos=new_pos(x,y,1,self.orientation)
            c=self.what_cell(npos)
            #if  self.near_margin():
             #   self.orientation=(self.orientation+180)%360          
              #  npos=new_pos(x,y,1,self.orientation)
               # c=self.what_cell(npos)
                
            if c not in cells: cells.append(c)
            attr=self.cell_attractiveness(cells)
            attr[c]*=2
            attr={k:v/float(sum(attr.values())) for k,v in attr.items()}
            
            
            if sum(attr.values())==1:
                target=attr.items()[np.random.choice(range(len(attr)),p=attr.values())][0]
                if target!=c: 
                    self.orientation=angle_cells(self.position,target)
                    print "d"
            self.wiggle(d)
            
        else:
            return False
            
        
    def has_food(self,radius,cells=None):
        
        if Dung.Instances.values()==[]: return False
        
        dung_patches=[]
        for d in Dung.Instances.values():
                       
            if euclid_dist(x1=d.position[0],x2=self.position[0],y1=d.position[1],y2=self.position[1])<radius:
                if cells!=None:
                    if self.what_cell(d.position) in cells: dung_patches.append(d)
                else:
                    dung_patches.append(d)
                
        if dung_patches==[]:
            return False    
        else:
            return dung_patches[np.random.randint(len(dung_patches))]
                                    
    def use_dung(self,dung, removed_dung,revolved_soil):
        if dung.amount<=0: 
            self.grid.remove_plume(dung.position,dung.reach)
            Dung.Instances.pop(dung.id)
        else:
            Dung.Functions[dung.id]['pedoturbation']+=revolved_soil
            Dung.Functions[dung.id]['incorporated_dung']+=removed_dung            
            dung.amount-=removed_dung
            #dung.incorporated_dung+= removed_dung
            #dung.pedoturbation+=revolved_soil
            self.energy+=removed_dung*8
            
            #self.rest_time.append(range(self.world.step,self.world.step+50))
                
        
        
        
    def action2(self,p={'rest':0.5,'move':0.5}):
        
        if self.world.step in self.rest_time or self.distance_covered>=self.max_dist: p={'rest':1.0,'move':0.0}
        if self.age<self.min_age: p={'rest':1.0,'move':0.0}
        if self.world.day not in self.active_days: p={'rest':1.0,'move':0.0}
        action=np.random.choice(p.keys(),p=p.values())
        
        food=self.has_food(self.perception_radius)
        #print action
        if action=='rest':
            self.energy-=self.energy_par['rest']
        elif food!=False and self.energy<self.energy_par['max_energy']:
             self.use_dung(food,self.functions()['removed_dung'],self.functions()['revolved_soil'])   

        else:            
            distance_to_move=np.random.normal(self.dist_par['mean'],self.dist_par['sd'],1)
            moore_neighborhood=self.grid.hood(self.what_cell(self.position),remove_center=False)
            possible_cells=self.cells_in_radius(cells=moore_neighborhood,radius=distance_to_move)
            cell_attr=self.cell_attractiveness(moore_neighborhood)
            
            #cells have different attractiveness
            if cell_attr.values().count(max(cell_attr.values()))<len(cell_attr):
                #move towards (one of) most attractive cells
                hv_cell=self.cell_with_higher_value(cell_attr)
                self.move_towards(d=distance_to_move,aim_position=hv_cell[0])
                #self.position=self.cell_center(hv_cell[0])
                #print "directional"
            else:
                #if they are all equally attractive, just wiggle    
                #print "random"
                self.wiggle(d=distance_to_move)
                
                
            self.energy-=self.energy_par['move']
            self.distance_covered+=distance_to_move
        self.track.append(self.position)
       
        
        
    
    def action(self,p={'rest':0.5,'move':0.5}):
        
        #if self.world.step in self.rest_time or self.distance_covered>=self.max_dist: p={'rest':1.0,'move':0.0}
        if self.age<self.min_age or self.world.day not in self.active_days:
            p={'rest':1.0,'move':0.0}
            #if self.world.day not in self.active_days: p={'rest':1.0,'move':0.0}
        else:
            food=self.has_food(self.perception_radius)
            if food!=False and self.energy<self.energy_par['max_energy']:
                 self.use_dung(food,self.functions()['removed_dung'],self.functions()['revolved_soil'])   
        #print action
        action=np.random.choice(p.keys(),p=p.values())
        if action=='rest':
            self.energy-=self.energy_par['rest']
               
        else:            
            distance_to_move=np.random.normal(self.dist_par['mean'],self.dist_par['sd'],1)
            moore_neighborhood=self.grid.hood(self.what_cell(self.position),remove_center=False)
            possible_cells=self.grid.circle(center=self.what_cell(self.position),radius=distance_to_move)
            suitable_habitats=self.habitat_suitability(possible_cells)
            cell_attr=self.cell_attractiveness(possible_cells)
            
            #cells have different attractiveness
            if cell_attr and cell_attr.values().count(max(cell_attr.values()))<len(cell_attr):
                #move towards (one of) most attractive cells
                hv_cell=self.cell_with_higher_value(cell_attr)
                #self.move_towards(d=distance_to_move,aim_position=hv_cell[0])
                self.position=self.cell_center(hv_cell[0])
                #print "directional"
            else:
                #if they are all equally attractive, just wiggle    
                #print "random"
                self.wiggle(d=distance_to_move,cells=suitable_habitats)
                
            if self.energy>20 and self.age>70 and self.world.day in self.breeding_days:
                partners=self.suitable_partners(suitable_habitats)
                if len(partners)>0:
                    self.breed(partner=np.random.choice(partners),cells=suitable_habitats)
                    
                
            self.energy-=self.energy_par['move']
            self.distance_covered+=distance_to_move
        self.track.append(self.position)
        
        
            
    def suitable_partners(self,cells):
            beetles=Beetle.BeetlesInArea(cells)
            partners=[]
            partners=[b for b in beetles if b.id!=self.id and b.sex!=self.sex and b.energy>20 and b.age>70]
            return partners
        
    def my_movement(self,d):
        self.wiggle(0)
        npos=new_pos(self.position[0],self.position[1],d,self.orientation)
        if not self.grid.in_bounds(npos): npos=self.grid.adjust_to_bounds(npos)
        if self.cell_suitability(npos):
            move(distance=d,angle=self.orientation)
        else:
            self.orientation=(self.orientation+180)%360
        
        
        
            
    def wiggle(self,d,sd=40,cells=None):
        if cells==None:cells=[]
        self.orientation+=np.random.randint(sd)
        self.orientation-=np.random.randint(sd)
        self.orientation=self.orientation%360
        npos=new_pos(self.position[0],self.position[1],d,self.orientation)
        if len(cells)>0 and npos not in cells:
            self.position=cells[np.random.randint(len(cells))]
        else:
             self.move(distance=d,angle=self.orientation)
            

    def functions(self):
        rd=0.5
        rs=0.3
        return {'removed_dung': rd, 'revolved_soil':rs}
    
    def near_margin(self):
        x=self.position[0]==0 or self.position[0]==self.grid.x_max
        y=self.position[1]==0 or self.position[1]==self.grid.y_max
        
        return x or y
    
        
    def cell_center(self,cell):
        return (floor(cell[0])+0.5, floor(cell[1])+0.5)
    
    def cells_in_radius(self,cells,radius, step=15):
        
        pos=self.position
        possible_cells=[]
        for alpha in xrange(0,360,step):
            p=self.what_cell((pos[0]+cos(radians(alpha))*radius,pos[1]+sin(radians(alpha))*radius))
            if p not in possible_cells: possible_cells.append(p) 
        
            
        return [cell for cell in possible_cells if cell in cells] 
        
    def habitat_suitability(self,cells):
        suitable=[]
        habitat=self.values_cells(dim="Habitat",cells=cells)    
        
        if self.energy<self.energy_par['hungry']:
            suitable=cells
        else:
            suitable=[c for c in cells if habitat[c]==max(habitat.values())]
        
        return suitable
            
        
        
        
        
    def cell_attractiveness(self,cells):
        habitat=self.values_cells(dim="Habitat",cells=cells)
        
        cues=self.values_cells(dim="Cues",cells=cells)
        sum_cues=float(sum(cues.values()))
        if sum_cues==0:
            cues={k:1./len(cells) for k,h in cues.items()}
        else:
            cues={k:v/sum_cues for k,v in cues.items()}
        
        
        if self.energy<self.energy_par['hungry']:
            habitat={k:1./len(cells) for k,h in habitat.items()}
            attractiveness={c:habitat[c]+(cues[c]) for c in cells}
            
                        
        elif self.energy_par['hungry']<self.energy<(self.energy_par['max_energy']-5):
            habitat={k:self.habitat_probs[h] for k,h in habitat.items()}
            sum_habitats=float(sum(habitat.values()))
            if sum_habitats==0:
                habitat={k:1./len(cells) for k,h in habitat.items()}
            else:    
                habitat={k:v/sum_habitats for k,v in habitat.items()}
            
            attractiveness={c:habitat[c]+(cues[c]) for c in cells}
        
        else:
            habitat={k:self.habitat_probs[h] for k,h in habitat.items()}
            sum_habitats=float(sum(habitat.values()))
            if sum_habitats==0:
                habitat={k:1./len(cells) for k,h in habitat.items()}
            else:    
                habitat={k:v/sum_habitats for k,v in habitat.items()}
            
            attractiveness={c:habitat[c] for c in cells}
            
            
        attractiveness={k:v/float(sum(attractiveness.values())) for k,v in attractiveness.items()}
        return attractiveness
        
    def values_cells(self,dim,cells):
        return {k:self.grid.value(dim=dim,cell=k) for k in cells}             
                    
    
    def values_moore(self,dim):
        h=self.grid.hood(cell=self.what_cell(self.position),remove_center=False)
        return {k:self.grid.value(dim=dim,cell=k) for k in h}
        
        
