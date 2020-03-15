import matplotlib.pyplot as plt
import numpy as np

class cluster:
    def __init__(self, step, hr=1.5, PowerW=100, R=500, f=1800,ht=32):
        self.resolution = step
        self.step = step
        self.hr = hr
        self.ht = ht
        self.ahre = 3.2*(np.log10(11.75*hr))**2 - 4.97
        self.len_grid=5*R
        self.subR = R*np.sqrt(3)/2
        self.power = 10*np.log10(PowerW/(10**-3))
        self.define_grid_axis(step,self.len_grid)
        self.R=R
        self.f=f

    def define_grid_axis(self,step,len_grid):
        axis = np.arange(step,len_grid-step,step)
        self.axis_x,self.axis_y=np.meshgrid(axis,axis)
    
    def get_cluster_info(self):
        return self.__dict__

    def create_cluster(self,X,Y,grid_len,subR,R):
        position1 = (X + 1j*Y) - (grid_len/2 + (grid_len/2)*1j)
        position2 = (X + 1j*Y) - (grid_len/2 + 2*subR*1j + (grid_len/2)*1j) 
        position3 = (X + 1j*Y) - (grid_len/2 - 2*subR*1j + (grid_len/2)*1j)
        position4 = (X + 1j*Y) - (grid_len/2 + (R + subR/2) + (subR*1j) + (grid_len/2)*1j)
        position5 = (X + 1j*Y) - (grid_len/2 + (R + subR/2) - (subR*1j) + (grid_len/2)*1j)
        position6 = (X + 1j*Y) - (grid_len/2 - (R + subR/2) + (subR*1j) + (grid_len/2)*1j) 
        position7 = (X + 1j*Y) - (grid_len/2 - (R + subR/2) - (subR*1j) + (grid_len/2)*1j) 
        final_cluster_grid_erbs_position=dict(erb1=position1,erb2=position2,erb3=position3,\
            erb4=position4,erb5=position5,erb6=position6,erb7=position7,)
        return final_cluster_grid_erbs_position
    
    def path_loss(self,data,pt,f,ht,ahre):
        result=[]
        for erb in range(1,len(data)+1):
            pDbm=pt - (46.3 + 33.9*np.log10(f) - 13.82*np.log10(ht) \
                - ahre + (44.9 - 6.55*np.log10(ht))*np.log10(np.absolute(data['erb%s'%(erb)]/1000)) + 3) 
            result.append(pDbm)            
        return result

    def union_erbs(self,result):
        new_color_map=np.zeros(result[0].shape) - 200
        for i  in range(0,len(result)):
            # print('result[%i],result[%i]'%(i,i+1))
            new_color_map=np.maximum(new_color_map,result[i])
        self.show_plot_color_map(new_color_map)
    
    def show_plot_color_map(self,color_map):
        print(color_map.shape)
        plt.pcolor(color_map,cmap='jet') 
        plt.grid()
        plt.colorbar(label="Potência em DB")
        plt.title(["Path Loss "])
        plt.xlabel("Distância em x10² M")
        plt.ylabel("Distância em x10² M")
        plt.show()


if __name__ == "__main__":
    Cluster=cluster(50)
    cluster_info=Cluster.get_cluster_info()
    data=Cluster.create_cluster(cluster_info['axis_x'],cluster_info['axis_y'],\
        cluster_info['len_grid'],cluster_info['subR'],cluster_info['R'])

    Path_loss=Cluster.path_loss(data,cluster_info['power'],cluster_info['f']\
        ,cluster_info['ht'],cluster_info['ahre'])

    Cluster.union_erbs(Path_loss)
 