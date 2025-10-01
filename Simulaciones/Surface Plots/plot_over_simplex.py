import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.tri as mtri
import matplotlib.animation as animation
import os


def coordXY(points): #Dados puntos en el símplex, devuelve las coordenadas XY en el plano.
    row, col= np.shape(points)
    coordxy=np.zeros((row, 2))
    for i in range(row):
        coordxy[i,:]=np.array([points[i, 1]+0.5*points[i,2],(np.sqrt(3)/2)*points[i, 2]])
    return coordxy


def surface_over_simplex(Xs,Ys,ZZs, col, col_borde,alfa,alfa_borde, elev=30, azim=-45):
    
    x=np.array([0,1,0.5])
    y=np.array([0,0,np.sqrt(3)/2])
    z=np.array([0,0,0])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev, azim)

    verts = [list(zip(x, y, z))]
     #grafico grilla en el piso
    ax.scatter3D(Xs, Ys, np.zeros(len(Xs)), marker='.', color="lightgrey", alpha=0.2, s=3)

    #Grafico el triángulo en el piso
    ax.add_collection3d(Poly3DCollection(verts, facecolor='whitesmoke', alpha=0.2))
    ax.set_xlim(0,1)
    ax.set_ylim(0,np.sqrt(3)/2)
    ax.set_zlim(0,1)
    ax.tick_params(axis='x', colors='grey') #color eje z
    #ax.set_zticks(np.linspace(0.0, 1.0, 11))
    ax.zaxis.set_tick_params(labelsize=6)
    x=np.array([0,1,0.5,0])
    y=np.array([0,0,np.sqrt(3)/2,0])
    z=np.array([0,0,0,0])
    ax.plot3D(x,y,z, color='grey', linewidth=0.75)
    #Para mover el eje z al otro lado
    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                        tmp_planes[0], tmp_planes[1], 
                        tmp_planes[4], tmp_planes[5])
    #Para ocultar los ejes x e y
    ax.grid(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().line.set_linewidth(0)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().line.set_linewidth(0)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide YZ Plane
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide XZ Plane
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    #Labels a los vértices
    ax.text(0.03, -0.1, 0, "$x_1$", color='black', size="small")
    ax.text(1.025, -0.05, 0, "$x_2$", color='black', size="small")
    ax.text(0.5, np.sqrt(3)/2+0.05, 0, "$x_3$", color='black', size="small")
   
    triang = mtri.Triangulation(Xs, Ys)
    #Saco los triángulos que quedan feos (depende el gráfico)
    posx=np.full(2500, fill_value="None")
    posx[np.arange(5*0,500*5+(5*0), 50)]=Xs[np.arange(5*0,500*5+(5*0), 50)]
    isBad = np.where(posx!="None", True, False)
    mask = np.all(isBad[triang.triangles],axis=1)
    triang.set_mask(mask)
    ############################################################
    #Hago el wireframe a mano
    for i in range(len(ZZs)):
        surf=ax.plot_trisurf(triang, ZZs[i], color=col[i], alpha=alfa[i] , linewidth=0.000001 , shade=False)
        for j in range(10):
            xl=Xs[np.arange(50*5*j,50*(5*j+1), 1)]
            yl=Ys[np.arange(50*5*j,50*(5*j+1), 1)]
            zl=(ZZs[i])[np.arange(50*5*j,50*(5*j+1), 1)]
            ax.plot3D(xl,yl,zl, color=col_borde[i],alpha=alfa_borde[i] ,linewidth=0.75)
        xl=Xs[np.arange(2450,2500, 1)]
        yl=Ys[np.arange(2450,2500, 1)]
        zl=(ZZs[i])[np.arange(2450,2500, 1)]
        ax.plot3D(xl,yl,zl, color=col_borde[i], alpha=alfa_borde[i],linewidth=0.75)
        for j in range(10):
            xl=Xs[np.arange(5*j,500*5+(5*j), 50)]
            yl=Ys[np.arange(5*j,500*5+(5*j), 50)]
            zl=(ZZs[i])[np.arange(5*j,500*5+(5*j), 50)]
            ax.plot3D(xl,yl,zl, color=col_borde[i],  alpha=alfa_borde[i],linewidth=0.75)
        xl=Xs[np.arange(49,2500, 50)]
        yl=Ys[np.arange(49,2500, 50)]
        zl=(ZZs[i])[np.arange(49,2500, 50)]
        ax.plot3D(xl,yl,zl, color=col_borde[i], alpha=alfa_borde[i], linewidth=0.75) 
    #fig.tight_layout()
    #plt.show()
    return fig


def surface_over_simplex_app(Xs,Ys,ZZs, elev=30, azim=-45):
    #Tiene el heatmap
    x=np.array([0,1,0.5])
    y=np.array([0,0,np.sqrt(3)/2])
    z=np.array([0,0,0])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev, azim)

    verts = [list(zip(x, y, z))]
     #grafico grilla en el piso
    ax.scatter3D(Xs, Ys, np.zeros(len(Xs)), marker='.', color="lightgrey", alpha=0.2, s=3)

    #Grafico el triángulo en el piso
    ax.add_collection3d(Poly3DCollection(verts, facecolor='whitesmoke', alpha=0.2))
    ax.set_xlim(0,1)
    ax.set_ylim(0,np.sqrt(3)/2)
    ax.set_zlim(0,1)
    ax.tick_params(axis='x', colors='grey') #color eje z
    #ax.set_zticks(np.linspace(0.0, 1.0, 11))
    ax.zaxis.set_tick_params(labelsize=6)
    x=np.array([0,1,0.5,0])
    y=np.array([0,0,np.sqrt(3)/2,0])
    z=np.array([0,0,0,0])
    ax.plot3D(x,y,z, color='grey', linewidth=0.75)
    #Para mover el eje z al otro lado
    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                        tmp_planes[0], tmp_planes[1], 
                        tmp_planes[4], tmp_planes[5])
    #Para ocultar los ejes x e y
    ax.grid(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().line.set_linewidth(0)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().line.set_linewidth(0)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide YZ Plane
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide XZ Plane
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    #Labels a los vértices
    ax.text(0.03, -0.1, 0, "$P$", color='black', size="small")
    ax.text(1.025, -0.05, 0, "$C$", color='black', size="small")
    ax.text(0.5, np.sqrt(3)/2+0.05, 0, "$L$", color='black', size="small")
   
    triang = mtri.Triangulation(Xs, Ys)
    #Saco los triángulos que quedan feos (depende el gráfico)
    posx=np.full(2500, fill_value="None")
    posx[np.arange(5*0,500*5+(5*0), 50)]=Xs[np.arange(5*0,500*5+(5*0), 50)]
    isBad = np.where(posx!="None", True, False)
    mask = np.all(isBad[triang.triangles],axis=1)
    triang.set_mask(mask)
    ############################################################
    #Hago el wireframe a mano
    for i in range(len(ZZs)):
        surf=ax.plot_trisurf(triang, ZZs[i], cmap='viridis' , linewidth=0.000001 , shade=False) 
    fig.tight_layout()
    fig.colorbar(surf, shrink=0.75)
    #plt.show()
    return fig

def surface_over_simplex_anim(Xs,Ys,ZZs, col, col_borde,alfa,alfa_borde, carpeta, archivo, elev=30):
    
    x=np.array([0,1,0.5])
    y=np.array([0,0,np.sqrt(3)/2])
    z=np.array([0,0,0])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    

    verts = [list(zip(x, y, z))]
     #grafico grilla en el piso
    ax.scatter3D(Xs, Ys, np.zeros(len(Xs)), marker='.', color="lightgrey", alpha=0.2, s=3)

    #Grafico el triángulo en el piso
    ax.add_collection3d(Poly3DCollection(verts, facecolor='whitesmoke', alpha=0.2))
    ax.set_xlim(0,1)
    ax.set_ylim(0,np.sqrt(3)/2)
    ax.set_zlim(0,1)
    ax.tick_params(axis='x', colors='grey') #color eje z
    #ax.set_zticks(np.linspace(0.0, 1.0, 11))
    ax.zaxis.set_tick_params(labelsize=6)
    x=np.array([0,1,0.5,0])
    y=np.array([0,0,np.sqrt(3)/2,0])
    z=np.array([0,0,0,0])
    ax.plot3D(x,y,z, color='grey', linewidth=0.75)
    #Para mover el eje z al otro lado
    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                        tmp_planes[0], tmp_planes[1], 
                        tmp_planes[4], tmp_planes[5])
    #Para ocultar los ejes x e y
    ax.grid(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().line.set_linewidth(0)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().line.set_linewidth(0)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide YZ Plane
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide XZ Plane
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    #Labels a los vértices
    ax.text(-0.05, -0.05, 0, "$x_1$", color='black', size="small")
    ax.text(1.025, -0.05, 0, "$x_2$", color='black', size="small")
    ax.text(0.5, np.sqrt(3)/2+0.05, 0, "$x_3$", color='black', size="small")
   
    triang = mtri.Triangulation(Xs, Ys)
    #Saco los triángulos que quedan feos (depende el gráfico)
    posx=np.full(2500, fill_value="None")
    posx[np.arange(5*0,500*5+(5*0), 50)]=Xs[np.arange(5*0,500*5+(5*0), 50)]
    isBad = np.where(posx!="None", True, False)
    mask = np.all(isBad[triang.triangles],axis=1)
    triang.set_mask(mask)
    ############################################################
    #Hago el wireframe a mano y grafico superficie
    verts = [list(zip(x, y, z))]
    #grafico grilla en el piso
    ax.scatter3D(Xs, Ys, np.zeros(len(Xs)), marker='.', color="lightgrey", alpha=0.2, s=3)

    #Grafico el triángulo en el piso
    ax.add_collection3d(Poly3DCollection(verts, facecolor='whitesmoke', alpha=0.2))
    ax.set_xlim(0,1)
    ax.set_ylim(0,np.sqrt(3)/2)
    ax.set_zlim(0,1)
    ax.tick_params(axis='x', colors='grey') #color eje z
    #ax.set_zticks(np.linspace(0.0, 1.0, 11))
    ax.zaxis.set_tick_params(labelsize=6)
    x=np.array([0,1,0.5,0])
    y=np.array([0,0,np.sqrt(3)/2,0])
    z=np.array([0,0,0,0])
    ax.plot3D(x,y,z, color='grey', linewidth=1)
    #Para mover el eje z al otro lado
    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                        tmp_planes[0], tmp_planes[1], 
                        tmp_planes[4], tmp_planes[5])
    #Para ocultar los ejes x e y
    ax.grid(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().line.set_linewidth(0)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().line.set_linewidth(0)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide YZ Plane
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide XZ Plane
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    #Labels a los vértices
    ax.text(-0.05, -0.05, 0, "$x_1$", color='black', size="small")
    ax.text(1.025, -0.05, 0, "$x_2$", color='black', size="small")
    ax.text(0.5, np.sqrt(3)/2+0.05, 0, "$x_3$", color='black', size="small")
    for i in range(len(ZZs)):
        ax.plot_trisurf(triang, ZZs[i], color=col[i], alpha=alfa[i] , linewidth=0.000001 , shade=False)
        for j in range(10):
            xl=Xs[np.arange(50*5*j,50*(5*j+1), 1)]
            yl=Ys[np.arange(50*5*j,50*(5*j+1), 1)]
            zl=(ZZs[i])[np.arange(50*5*j,50*(5*j+1), 1)]
            ax.plot3D(xl,yl,zl, color=col_borde[i],alpha=alfa_borde[i] ,linewidth=0.75)
        xl=Xs[np.arange(2450,2500, 1)]
        yl=Ys[np.arange(2450,2500, 1)]
        zl=(ZZs[i])[np.arange(2450,2500, 1)]
        ax.plot3D(xl,yl,zl, color=col_borde[i], alpha=alfa_borde[i],linewidth=0.75)
        for j in range(10):
            xl=Xs[np.arange(5*j,500*5+(5*j), 50)]
            yl=Ys[np.arange(5*j,500*5+(5*j), 50)]
            zl=(ZZs[i])[np.arange(5*j,500*5+(5*j), 50)]
            ax.plot3D(xl,yl,zl, color=col_borde[i],  alpha=alfa_borde[i],linewidth=0.75)
        xl=Xs[np.arange(49,2500, 50)]
        yl=Ys[np.arange(49,2500, 50)]
        zl=(ZZs[i])[np.arange(49,2500, 50)]
        ax.plot3D(xl,yl,zl, color=col_borde[i], alpha=alfa_borde[i], linewidth=0.75) 
    def animate(n):
        ax.view_init(elev, n)
        ax.set_title(str("Ángulo Azimut = ")+str(n))
        
        return fig
    os.chdir(carpeta)

    anim = animation.FuncAnimation(fig = fig, func = animate, frames = 360, interval=20, repeat = False)
    
    writergif = animation.PillowWriter(fps=20)
    anim.save(carpeta+str("\\")+str(archivo)+str(".gif"),writer=writergif)
    


def surface_over_simplex_heatmap(Xs,Ys,Zs, elev=30, azim=-45):
    
    x=np.array([0,1,0.5])
    y=np.array([0,0,np.sqrt(3)/2])
    z=np.array([0,0,0])


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev, azim)

    verts = [list(zip(x, y, z))]
    #grafico grilla en el piso
    ax.scatter3D(Xs, Ys, np.zeros(len(Xs)), marker='o', color="grey", alpha=0.3, s=2)

    #Grafico el triángulo en el piso
    ax.add_collection3d(Poly3DCollection(verts, facecolor='whitesmoke', alpha=0.2))
    ax.set_xlim(0,1)
    ax.set_ylim(0,np.sqrt(3)/2)
    ax.set_zlim(0,1)
    #ax.set_zticks(np.linspace(0.0, 1.0, 11))
    ax.zaxis.set_tick_params(labelsize=6)
    x=np.array([0,1,0.5,0])
    y=np.array([0,0,np.sqrt(3)/2,0])
    z=np.array([0,0,0,0])
    ax.plot3D(x,y,z, color='grey', linewidth=0.75)
    #Para mover el eje z al otro lado
    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                        tmp_planes[0], tmp_planes[1], 
                        tmp_planes[4], tmp_planes[5])
    #Para ocultar los ejes x e y
    ax.grid(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().line.set_linewidth(0)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().line.set_linewidth(0)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide YZ Plane
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) # Hide XZ Plane
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    #Labels a los vértices
    ax.text(0.03, -0.1, 0, "$x_1$", color='black')
    ax.text(1.025, -0.025, 0, "$x_2$", color='black')
    ax.text(0.5, np.sqrt(3)/2+0.05, 0, "$x_3$", color='black')
    
    triang = mtri.Triangulation(Xs, Ys)
    #Saco los triángulos que quedan feos (depende el gráfico)############
    posx=np.full(2500, fill_value="None")
    posx[np.arange(5*0,500*5+(5*0), 50)]=Xs[np.arange(5*0,500*5+(5*0), 50)]
    isBad = np.where(posx!="None", True, False)
    mask = np.all(isBad[triang.triangles],axis=1)
    triang.set_mask(mask)
    ##################################################################
    surf=ax.plot_trisurf(triang, Zs, cmap='viridis', alpha=0.95, shade=False, linewidth=0.00001)
    fig.tight_layout()
    fig.colorbar(surf, shrink=0.75)
    
    return fig



