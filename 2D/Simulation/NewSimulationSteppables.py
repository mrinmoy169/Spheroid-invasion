from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import random
import numpy as np
from collections import Counter

class NewSimulationSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        self.force_mod = -40.01  #active motility force (value of mu)

    def start(self):
        #creating the initial setup for splenocytes
        for ix in np.arange(0,self.dim.x,15):
            for iy in np.arange(0,self.dim.y,15):
                if np.sqrt((ix-self.dim.x/2)**2+(iy-self.dim.y/2)**2)>105:
                    self.cell_field[int(ix):int(ix)+5, int(iy):int(iy)+5, 0] = self.new_cell(self.T)
        
        for cell in self.cell_list:
            if cell.type == self.E:
                cell.targetVolume = 25
                cell.lambdaVolume = 1.0
                cell.targetSurface = 25
                cell.lambdaSurface = 50.0
            if cell.type == self.T:
                cell.targetVolume = 25
                cell.lambdaVolume = 1.0
                cell.targetSurface = 25
                cell.lambdaSurface = 1.0
            
                angle = np.random.normal(0,2*np.pi)
                cell.lambdaVecX= self.force_mod*np.cos(angle)
                cell.lambdaVecY= self.force_mod*np.sin(angle)

            cell.dict["postauX"]=cell.xCOM
            cell.dict["postauY"]=cell.yCOM
            cell.dict["vx"]=0.0
            cell.dict["vy"]=0.0

        global no_of_T
        no_of_T=[]
    def step(self,mcs):  
        for cell in self.cell_list_by_type(self.T):
            
            if mcs%10==0: #10 is the persistence time of cell polarity (tau)
                dx=(cell.xCOM - cell.dict["postauX"])
                dy=(cell.yCOM - cell.dict["postauY"])
                #pbc correction
                if dx>0.5*self.dim.x:
                    dx-=self.dim.x
                if dx<=-0.5*self.dim.x:
                    dx+=self.dim.x
                if dy>0.5*self.dim.y:
                    dy-=self.dim.y
                if dy<=-0.5*self.dim.y:
                    dy+=self.dim.y
                vx = dx/10
                vy = dy/10
                cell.dict["postauX"]=cell.xCOM
                cell.dict["postauY"]=cell.yCOM
                cell.dict["vx"]=vx
                cell.dict["vy"]=vy
                norm_v = np.linalg.norm([vx,vy])
                if norm_v != 0:
                    vx /= norm_v
                    vy /= norm_v
                
                cell.lambdaVecX = self.force_mod*vx
                cell.lambdaVecY = self.force_mod*vy
        
        for cell in self.cell_list_by_type(self.T):
            area=0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor is not None and neighbor.type == self.E:
                    area += common_surface_area
            if area > 0.95*cell.surface: #0.95 is the contact dependent growth parameter f
                cell.targetVolume += 25/20000 #20000 is the doubling time of splenocytes
                cell.targetSurface = 5*np.sqrt(cell.targetVolume)
  
        if mcs%100==0:
            cluster1=[]
            cluster2=[]
            cluster1_cell=[]
            cluster2_cell=[]
            for cell in self.cell_list_by_type(self.T):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor is not None and neighbor.type == self.E:
                        cluster1.append(cell.id)
                        cluster1_cell.append(cell)
                    if neighbor is None:
                        cluster2.append(cell.id)
                        cluster2_cell.append(cell)
            cluster1_unique=np.unique(cluster1)
            cluster1_unique_cell=list(Counter(cluster1_cell))
            cluster3=[x for x in cluster2 if x not in cluster1_unique]
            cluster3_cell=[x for x in cluster2_cell if x not in cluster1_unique_cell]
            Tin = len(cluster1_unique)
            Tout = len(cluster3)
            Tin_area = 0
            for cell in cluster1_unique_cell:
                Tin_area += cell.volume
            Tout_area = 0
            for cell in cluster3_cell:
                Tout_area += cell.volume
            no_of_T.append((mcs,Tin,Tout,Tin_area,Tout_area))
            
        output_path = Path(self.output_dir).joinpath('no_of_T.dat')
        np.savetxt(output_path,no_of_T,fmt='%12.2f') #saving the splenocyte pixels inside and outside the tumor
            
#calculating the splenocite positions
        if mcs%1000==0:
            list_r=[]
            list_r1=[]
            for cell in self.cell_list_by_type(self.T):
                list_r.append(np.sqrt((cell.xCOM-self.dim.x/2)**2+(cell.yCOM-self.dim.y/2)**2))
                pixel_list = self.get_cell_pixel_list(cell)
                for pixel_tracker_data in pixel_list:
                    pt=pixel_tracker_data.pixel
                    list_r1.append(np.sqrt((pt.x-self.dim.x/2)**2+(pt.y-self.dim.y/2)**2))

#calculating the boundary of the tumor
        if mcs%1000==0:
            list_T1=[]
            list_T2=[]
            for cell in self.cell_list_by_type(self.T):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor is not None and neighbor.type == self.E:
                        list_T1.append(cell)   
                    if neighbor is None:
                        list_T2.append(cell)
            list_T1_unique=list(Counter(list_T1))
            list_T3=[x for x in list_T1_unique if x not in list_T2]
            list_T4=[x for x in list_T2 if x not in list_T1_unique]
            list_T_cell=[x for x in list_T1 if x not in list_T3]
            list_T=[]
            for cell in list_T_cell:
                list_T.append(cell.id)   
            list_E1=[]
            list_E2=[]
            for cell in self.cell_list_by_type(self.E):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor is not None and neighbor.id in list_T:
                        list_E1.append(cell)
                    if neighbor is None:
                        list_E2.append(cell)
            list_E3=list_E1+list_E2
            list_E=list(Counter(list_E3))
  
                
        if mcs%1000==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('boundary' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for item in list_E:
                        fout.write('{} {}\n'.format(item.xCOM, item.yCOM))
        
        if mcs%1000==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('step_' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for cell in self.cell_list:
                        fout.write('{} {} {} {} {} {}\n'.format(cell.id, cell.type, cell.xCOM, cell.yCOM, cell.volume, cell.surface))

        if mcs%1000==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('T_r_' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for item in list_r:
                        fout.write('{}\n'.format(item))
                        
        if mcs%1000==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('T_r1_' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for item in list_r1:
                        fout.write('{}\n'.format(item))

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def start(self):
        
        self.initial_volume = 25
            
    def step(self, mcs):

        cells_to_divide = [] # a list of cells that have to go through division
        for cell in self.cell_list_by_type(self.T):
            if cell.volume >= 2*self.initial_volume:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0
        self.parent_cell.targetSurface = 5*np.sqrt(self.parent_cell.targetVolume)
        self.clone_parent_2_child()    
