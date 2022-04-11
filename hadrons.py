import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mpl3
import numpy as np
from tabulate import tabulate
import os
import sys

class quark():
    def __init__(self, name, title, Q=0, B=1/3, Iz=0, S=0, C=0, Bt=0, T=0):
        self.name = name
        self.title = title
        self.Q    = Q     # electric charge
        self.B    = B     # baryon number
        self.Iz = Iz    # 3rd component of Isospin
        self.S    = S     # strangeness
        self.C    = C     # charmness
        self.Bt = Bt    # bottomness
        self.T    = T     # topness
        self.Y    = self.B + self.S + self.C + self.Bt + self.T
        if abs(self.Y - 2*(self.Q-self.Iz))>1e-6:
            raise RuntimeError('Creation of quark '+name+' failed. Hypercharge inconsistent '+str(self.Y)+' '+str(2*(self.Q-self.Iz)))

    def anti(self,name='',title=''):
        if name=='': name = 'anti-'+self.name
        return quark(name, title, -self.Q, -self.B, -self.Iz, -self.S, -self.C, -self.Bt, -self.T)

    def __str__(self):
        return 'quark({:s},Iz={:3.1f},S={:d},C={:d},Bt={:d},T={:d})'.format(self.name, self.Iz, self.S, self.C, self.Bt, self.T)

class state():
    def __init__(self,name,title,quarks):
        self.name     = name
        self.title    = title
        self.quarks = quarks
        self.Q     = sum( [ q.Q    for q in quarks ] )
        self.B     = sum( [ q.B    for q in quarks ] )
        self.Iz    = sum( [ q.Iz for q in quarks ] )
        self.S     = sum( [ q.S    for q in quarks ] )
        self.C     = sum( [ q.C    for q in quarks ] )
        self.Bt    = sum( [ q.Bt for q in quarks ] )
        self.T     = sum( [ q.T    for q in quarks ] )
        self.Y    = self.B + self.S + self.C + self.Bt + self.T

        if abs(self.Y - 2*(self.Q-self.Iz))>1e-6:
            raise RuntimeError('Creation of quark '+name+' failed. Hypercharge inconsistent '+str(self.Y)+' '+str(2*(self.Q-self.Iz)))

    def anti(self,name='',title=''):
        return state(name,title, [q.anti() for q in self.quarks])

    def __str__(self):
        content = ','.join( [ q.name for q in self.quarks ] )
        return '{:s}({:s}) Iz={:3.1f},S={:d},C={:d},B={:d})'.format(self.name,content,self.Iz, self.S, self.C, self.B)

def make_anti(parts):
    return ( (-a,-b) for a,b in parts )

def meson(a,b):
    return (a[0]+b[0],a[1]+b[1])

def plot(states, ax=None, rad=0.15, hexc=(0,0), hexr=1, tri=None):

    ax = ax or plt.gca()

    # min, max
    ymin = min( [state.S for state in states] )
    ymax = max( [state.S for state in states] )

    xmin = min( [state.Iz for state in states] )
    xmax = max( [state.Iz for state in states] )

    # axis set up
    ax.set_xticks( np.linspace(xmin,xmax,5) )
    ax.tick_params(axis='x', labelsize=14)
    ax.set_xlabel('Isospin, $I_z$', fontsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_ylabel('Strangeness, $S$', fontsize=14)

    # buffer
    buff = 0.2
    xmin -= buff
    xmax += buff
    ymin -= buff
    ymax += buff

    # strange lines
    for s in np.linspace(-2,2,5):
        ax.plot( (xmin, xmax), (s,s), 'k--', lw=1 )

    # isospin lines
    #for i in np.linspace(-1,1,5):
        #ax.plot( (i,i), (ymin,ymax), 'k--', lw=1 )

    # electric charge lines
    for q in np.linspace(-2,2,5):
        ax.plot( (xmin+hexc[1]/2, q-ymin/2+hexc[1]/2), (2*q-2*xmin, ymin), 'k--', lw=1 )

    # draw the hexes
    if hexr is not None:
        hx = np.array([ -hexr, -hexr/2, hexr/2, hexr,    hexr/2, -hexr/2, -hexr ])
        hy = np.array([         0,    hexr    , hexr    ,        0, -hexr    , -hexr    ,         0 ])
        hx += hexc[0]
        hy += hexc[1]
        ax.plot(hx,hy,'r-', lw=3)

    # draw the triangle
    if tri is not None:
        ax.plot( tri[0], tri[1], 'r-', lw=3 )

    # draw the states
    coords = [ (state.Iz, state.S) for state in states ]

    # set circle radius
    rep = 0

    for i, state in enumerate(states):
        (x, y) = (state.Iz, state.S)
        # check how many entries for that coord
        n = coords.count( coords[i] )
        if n==2:
            if rep==0: x+=rad
            else: x-=rad
            rep+=1
        if n==3:
            if rep==0:
                y += 0.866*rad
            elif rep==1:
                x -= rad
                y -= 0.866*rad
            elif rep==2:
                x += rad
                y -= 0.866*rad
            rep+=1

        circ = plt.Circle( (x,y), 0.15, ec='r', fc='lightblue', lw=3, zorder=10 )
        ax.add_patch(circ)
        if state.title!='':
            ax.text( x, y, state.title, ha='center', va='center', fontsize=16, zorder=20)

    #ax.axis('off')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    #ax.set_aspect('equal')#,'box')

def hex3d( hexc=(0,0), hexr=1, z=0 ):

    hx = np.array([ -hexr, -hexr/2, hexr/2, hexr,    hexr/2, -hexr/2 ])
    hy = np.array([         0,    hexr    , hexr    ,        0, -hexr    , -hexr     ])
    hx += hexc[0]
    hy += hexc[1]
    hz = np.full_like(hx,z)

    vtxs = np.column_stack( (hx,hy,hz) )

    return mpl3.art3d.Poly3DCollection( [vtxs] )

def tri3d( v1=(-0.5,-1), v2=(0,1), v3=(0.5,-1), z=1 ):

    hx = np.array( [v1[0],v2[0],v3[0]] )
    hy = np.array( [v1[1],v2[1],v3[1]] )
    hz = np.full_like(hx,z)

    vtxs = np.column_stack( (hx,hy,hz) )

    return mpl3.art3d.Poly3DCollection( [vtxs] )

def plot3d(states=None, ax=None, mes=0, bar=None, content=True):

    ax = ax or plt.gca()

    if mes is not None:
        # hex
        shape = hex3d()
        shape.set_color('lightblue')
        ax.add_collection(shape)

        # tri
        tri1 = tri3d( (-0.5,0), (0, 1), (0.5,0), z= 1 )
        tri2 = tri3d( (-0.5,0), (0,-1), (0.5,0), z=-1 )
        tri1.set_color('lightblue')
        tri2.set_color('lightblue')
        ax.add_collection(tri1)
        ax.add_collection(tri2)

        # lines
        ax.plot( (-0.5,-1    ,-0.5), (0, 0, 0), (1,0,-1), 'k-', lw=1 )
        ax.plot( ( 0.5, 1    , 0.5), (0, 0, 0), (1,0,-1), 'k-', lw=1 )
        ax.plot( ( 0    ,-0.5,-0.5), (1, 1, 0), (1,0,-1), 'k-', lw=1 )
        ax.plot( ( 0    , 0.5, 0.5), (1, 1, 0), (1,0,-1), 'k-', lw=1 )
        ax.plot( (-0.5,-0.5, 0    ), (0,-1,-1), (1,0,-1), 'k-', lw=1, zorder=50 )
        ax.plot( ( 0.5, 0.5, 0    ), (0,-1,-1), (1,0,-1), 'k-', lw=1, zorder=50 )

    if bar is not None:
        if bar==0:
            # hex
            shape = hex3d( hexc=(0,-1), hexr=1, z=0 )
            shape.set_color('lightblue')
            ax.add_collection(shape)
            # tri
            tri1 = tri3d( (0,-2), (-    1,0), (    1,0), z= 1 )
            tri2 = tri3d( (0,-1), (-0.5,0), (0.5,0), z= 2 )
            tri1.set_color('lightblue')
            tri2.set_color('lightblue')
            ax.add_collection(tri1)
            ax.add_collection(tri2)
            # lines
            ax.plot( (-0.5, -1, -1), (0,0,-1), (2,1,0), 'k-', lw=1, zorder=10 )
            ax.plot( (-1,-0.5), (0,0), (1,0), 'k-', lw=1, zorder=10 )
            ax.plot( (0.5, 1, 1), (0,0,-1), (2,1,0), 'k-', lw=1, zorder=10 )
            ax.plot( (1,0.5), (0,0), (1,0), 'k-', lw=1, zorder=10 )
            ax.plot( (0,0,-0.5), (-1,-2,-2),(2,1,0),'k-',lw=1,zorder=10)
            ax.plot( (0,0.5),(-2,-2),(1,0),'k-',lw=1,zorder=10)


        if bar==1:
            # tri
            tri1 = tri3d( (0,-3), (-1.5,0), (1.5,0), z= 0 )
            tri2 = tri3d( (0,-2), (-1    ,0), (1    ,0), z= 1 )
            tri3 = tri3d( (0,-1), (-0.5,0), (0.5,0), z= 2 )
            tri1.set_color('lightblue')
            tri2.set_color('lightblue')
            tri3.set_color('lightblue')
            tri1.set_edgecolor('k')
            tri2.set_edgecolor('k')
            tri3.set_edgecolor('k')
            ax.add_collection(tri1)
            ax.add_collection(tri2)
            ax.add_collection(tri3)
            # lines
            ax.plot( (0, 0    ), (0,-3), (3,0), 'k-', lw=1, zorder=10 )
            ax.plot( (0,-1.5), (0, 0), (3,0), 'k-', lw=1, zorder=10 )
            ax.plot( (0, 1.5), (0, 0), (3,0), 'k-', lw=1, zorder=10 )

    # add states
    #x = [ state.Iz for state in states ]
    #y = [ state.S    for state in states ]
    #z = [ state.C    for state in states ]
    coords = [ (state.Iz, state.S, state.C) for state in states ]
    rad = 0.025
    ysc = 5
    rep = 0
    xps = []
    yps = []
    zps = []
    for i, state in enumerate(states):
        (x, y, z) = (state.Iz, state.S, state.C)
        # check how many entries for that coord
        n = coords.count( coords[i] )
        if n==2:
            if rep%2==0: x+=rad
            else: x-=rad
            rep+=1
        if n==3:
            if rep==0:
                y += 5*rad
            elif rep==1:
                x -= rad
                y -= 5*rad
            elif rep==2:
                x += rad
                y -= 5*rad
            rep+=1
        if n==4:
            if rep==0:
                x -= rad
                y += 5*rad
            elif rep==1:
                x += rad
                y += 5*rad
            elif rep==2:
                x -= rad
                y -= 5*rad
            elif rep==3:
                x += rad
                y -= 5*rad
            rep+=1

        xps.append( x )
        yps.append( y )
        zps.append( z )

    #ax.scatter( xps, yps, zps, marker='o', c='k', alpha=1, zorder=10 )
    #ax.scatter( xps, yps, zps, marker='o', c='k', alpha=1, zorder=20 )
    ax.plot( xps, yps, zps, 'ko', zorder=20)

    # state labels
    for i, state in enumerate(states):
        xsgn = 1
        zsgn = 1
        if mes is not None:
            if mes==0:
                xsgn = 0.7 if xps[i]>=0 else -2.3
                zsgn = -2 if yps[i]<0 or zps[i]<0 else 0.8
                if yps[i]==-1 and xps[i]==0.5: xsgn = -0.3
                if yps[i]<0 and yps[i]>-1: xsgn -= 0.6
                if yps[i]<0 and yps[i]>-1: zsgn += 0.4
                if yps[i]<0 and zps[i]<0 : zsgn += 0.4
                if xps[i]==0 and yps[i]==-1:
                    xsgn += 0.5
                    zsgn += 2.5
                if xps[i]==1: xsgn += 0.4
            elif mes==1:
                xsgn = 0.7 if xps[i]>=0 else -2.3
                zsgn = -2.2 if yps[i]<0 or zps[i]<0 else 0.8
                if yps[i]==-1 and xps[i]==0.5: xsgn = -0.3
                if yps[i]<0 and yps[i]>-1: xsgn -= 0.6
                if yps[i]<0 and yps[i]>-1: zsgn += 0.4
                if yps[i]<0 and zps[i]<0 : zsgn += 0.6
                if xps[i]==0.5: xsgn -= 0.3
                if xps[i]==0.5 and yps[i]==-1:
                    xsgn -= 0.2
                    zsgn -= 0.2
                if xps[i]==0 and yps[i]==-1:
                    xsgn += 0.5
                    zsgn += 2.5
                if xps[i]==1: xsgn += 0.4
        if bar is not None:
            if bar==0:
                xsgn = 0.8 if xps[i]>=0 else -3.0
                zsgn = 1 if yps[i]==0 else -2
                zsgn = 1.5 if yps[i]==0 and abs(xps[i])<1.5 else zsgn
                if yps[i]==0 and zps[i]==0:
                    xsgn *= 1.2
                    zsgn -= 1
                if yps[i]==0 and zps[i]==1 and xps[i]<0 and xps[i]>-1: xsgn += 0.4
                if yps[i]==-1 and zps[i]==1 and xps[i]<0 and xps[i]>-1: xsgn += 0.8
                if yps[i]==-2 and zps[i]==1:
                    xsgn += 0.6
                    zsgn += 0.6
                if yps[i]==-2 and zps[i]==0:
                    zsgn += 3
                    xsgn *= 0.5
                if yps[i]==-1 and zps[i]==0 and xps[i]>0: xsgn -= 0.4
            elif bar==1:
                xsgn = 1 if xps[i]>=0 else -5.5
                zsgn = 1 if yps[i]==0 else -2
                zsgn = 1.5 if yps[i]==0 and abs(xps[i])<1.5 else zsgn
                if xps[i]==0 and yps[i]==-3 and zps[i]==0: xsgn += 0.4
                if xps[i]==1.5 and yps[i]==0 and zps[i]==0: zsgn += 0.4

        ax.text( xps[i]+xsgn*0.05, yps[i]+0.05, zps[i]+zsgn*0.08, state.title, ha='left', va='center', zorder=20, fontsize=14 )
        if content:
            qstr = ''.join( [q.title for q in state.quarks ] )
            xp = xps[i]
            yp = yps[i]
            zp = zps[i]
            if mes is not None:
                if abs(xps[i])<0.5 and abs(yps[i])<1: continue
                if xp<0: xp += 0.05
                if xp>0: xp -= 0.12
                if zp==1 and abs(yp)==1:
                    xp -= 0.05
                    zp -= 0.1
                if zp==-1 and abs(yp)==1:
                    xp -= 0.03
                    zp += 0.1
            if bar is not None:
                if bar==0:
                    if abs(xps[i])<0.5 and abs(yps[i])==1 and zps[i]==0: continue
                    if xp<0: xp += 0.05
                    if xp>0: xp -= 0.16
                    if xp==0:
                        xp += 0.02
                        zp += 0.08
                    if yp==0 and zp==1 and abs(xp)<0.5:
                        zp -= 0.12
                        xp += 0.04
                    if yps[i]==-2 and zps[i]==0:
                        zp += 0.04
                    if yps[i]==0 and zps[i]==0:
                        zp += 0.04
                elif bar==1:
                    if xp<0: xp += 0.06
                    if xp>0: xp -= 0.23
                    if abs(xp)+abs(zp)>1: zp += 0.05
                    if yps[i]==0 and zps[i]==0 and abs(xps[i])==0.5: zp -= 0.12
                    if yps[i]==0 and zps[i]==1 and xps[i]==0: zp -= 0.12
                    if yps[i]==-1 and xps[i]==0 and zps[i]==0: zp += 0.1
                    if abs(yp)+abs(zp)>2 and xp==0:
                        xp += 0.02
                        zp += 0.08
                        if zp<2: zp += 0.04
                    if zp>=3:
                        xp -= 0.2
                        zp -= 0.04

            ax.text( xp, yp, zp, qstr, ha='left', va='center', zorder=20, fontsize=8 )

    # axes
    ax.set_xlabel('Isospin, $I_z$', fontsize=12)
    ax.set_ylabel('Strangeness, $S$', fontsize=12)
    ax.set_zlabel('Charmness, $D$', fontsize=12)
    if mes is not None:
        ax.margins(x=0,y=0,z=0)
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(-1,1)
        ax.set_xticks( np.linspace(-1,1,5) )
        ax.set_yticks( np.linspace(-1,1,3) )
        ax.set_zticks( np.linspace(-1,1,3) )
    if bar is not None:
        if bar == 0:
            ax.set_xlim(-1,1)
            ax.set_ylim(-2,0)
            ax.set_zlim(-0,2)
            ax.set_xticks( np.linspace(-1,1,5) )
            ax.set_yticks( np.linspace(-2,0,3) )
            ax.set_zticks( np.linspace( 0,2,3) )
        if bar == 1:
            ax.set_xlim(-1.5,1.5)
            ax.set_ylim(-3,0)
            ax.set_zlim(-0,3)
            ax.set_xticks( np.linspace(-1.5,1.5,7) )
            ax.set_yticks( np.linspace(-3,0,4) )
            ax.set_zticks( np.linspace( 0,3,4) )

    ax.view_init(elev=18,azim=-84)

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =    (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =    (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =    (1,1,1,0)
    # draw grid lines by hand because matplotlib adds a buffer
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()
    xlf    = xlim[0]
    xrt    = xlim[1]
    ybk    = ylim[1]
    yfr    = ylim[0]
    zpt    = zlim[0]
    for v in ax.get_xticks():
        ax.plot( (v,v), (ybk,ybk), zlim, lw=1, c='0.7', alpha=0.5 )
        ax.plot( (v,v), (yfr,yfr), zlim, lw=1, c='0.7', alpha=0.5, zorder=20 )
        for z in ax.get_zticks():
            ax.plot( (v,v), ylim, (z,z), lw=1, c='0.7', alpha=0.5 )
    for v in ax.get_yticks():
        for z in ax.get_zticks():
            ax.plot( xlim, (v,v), (z,z), lw=1, c='0.7', alpha=0.5 )
        if v==ax.get_yticks()[0] or v==ax.get_yticks()[-1]: continue
        ax.plot( (xlf,xlf), (v,v), zlim, lw=1, c='0.7', alpha=0.5 )
        ax.plot( (xrt,xrt), (v,v), zlim, lw=1, c='0.7', alpha=0.5 )
    for v in ax.get_zticks():
        ax.plot( xlim, (yfr,yfr), (v,v), lw=1, c='0.7', alpha=0.5, zorder=20 )

    # draw axis by hand because matplotlib adds buffer
    # x axis
    ax.plot( xlim, (yfr,yfr), (zpt,zpt), 'k-' , zorder=20)
    for v in ax.get_xticks():
        ax.plot( (v,v), (yfr,yfr), (zpt,zpt-0.02), 'k-' )
        ax.text( v, yfr, zpt-0.05*(zlim[1]-zlim[0]), str(v), ha='center', va='center' )
    ax.text( xlim[0] + (xlim[1]-xlim[0])/2, yfr, zpt-0.13*(zlim[1]-zlim[0]), 'Isospin, $I_z$', ha='center', va='center' )
    # y axis
    ax.plot( (xrt,xrt), ylim, (zpt,zpt), 'k-' )
    for v in ax.get_yticks():
        ax.plot( (xrt,xrt+0.02), (v,v), (zpt,zpt), 'k-' )
        ax.text( xrt+0.05*(xlim[1]-xlim[0]), v, zpt, str(v), ha='center', va='center' )
    ax.text( xrt+0.13*(xlim[1]-xlim[0]), ylim[0] + (ylim[1]-ylim[0])/2, zpt, 'Strangeness, $S$', ha='center', va='center', zdir='y' )
    # z axis
    ax.plot( (xrt,xrt), (ybk,ybk), zlim, 'k-' )
    for v in ax.get_zticks():
        ax.plot( (xrt,xrt+0.02), (ybk,ybk), (v,v), 'k-' )
        if v==ax.get_zticks()[0]: continue
        ax.text( xrt+0.05*(xlim[1]-xlim[0]), ybk, v, str(v), ha='center', va='center' )
    ax.text( xrt+0.13*(xlim[1]-xlim[0]), ybk, zlim[0] + (zlim[1]-zlim[0])/2, 'Charmness, $C$', ha='center', va='center', zdir='z' )


    ax.axis('off')

def print_quarks(quarks):
    rows    = []
    heads = ['Name','Q','B','Iz','Y','S','C','B\'','T']
    for key, item in quarks.items():
        Q = '{: d}/3'.format( round(3*item.Q) )
        B = '{: d}/3'.format( round(3*item.B) )
        Iz = str(item.Iz) if item.Iz==0 else '{: d}/2'.format( round(2*item.Iz) )
        Y = '{: d}/3'.format( round(3*item.Y) )
        rows.append( [ item.name, Q, B, Iz, Y, item.S, item.C, item.Bt, item.T ] )
    print('\033[1m Quarks \033[0m')
    print(tabulate(rows,tablefmt='pretty',headers=heads,colalign=('left','right','right','right','right','right','right','right','right')))

def print_states(states, head='States'):
    rows    = []
    heads = ['Name','Quarks','Q','B','Iz','Y','S','C','B\'','T']
    for key, item in states.items():
        quarks = '(' + ','.join( [qs.name for qs in item.quarks] ) + ')'
        Q = '{: d}/3'.format( round(3*item.Q) )
        B = '{: d}/3'.format( round(3*item.B) )
        Iz = str(item.Iz) if item.Iz==0 else '{: d}/2'.format( round(2*item.Iz) )
        Y = '{: d}/3'.format( round(3*item.Y) )
        rows.append( [ item.name, quarks, Q, B, Iz, Y, item.S, item.C, item.Bt, item.T ] )
    print('\033[1m %s \033[0m'%head)
    print(tabulate(rows,tablefmt='pretty',headers=heads,colalign=('left','left','right','right','right','right','right','right','right','right')))

if __name__ == "__main__":

    os.system('mkdir -p pdf')
    os.system('mkdir -p png')

    # make quarks
    u = quark( 'u', '$u$', Q =    2/3, Iz =    0.5 )
    d = quark( 'd', '$d$', Q = -1/3, Iz = -0.5 )
    s = quark( 's', '$s$', Q = -1/3, S    = -1     )
    c = quark( 'c', '$c$', Q =    2/3, C    =    1     )
    b = quark( 'b', '$b$', Q = -1/3, Bt = -1     )
    t = quark( 't', '$t$', Q =    2/3, T    =    1     )

    ubar = u.anti('ubar', r'$\bar{u}$')
    dbar = d.anti('dbar', r'$\bar{d}$')
    sbar = s.anti('sbar', r'$\bar{s}$')
    cbar = c.anti('cbar', r'$\bar{c}$')
    bbar = b.anti('bbar', r'$\bar{b}$')
    tbar = t.anti('tbar', r'$\bar{t}$')

    quarks = { 'u': u, 'd': d, 's': s, 'c': c, 'b': b, 't': t,
               'ubar': ubar, 'dbar': dbar, 'sbar': sbar, 'cbar': cbar, 'bbar': bbar, 'tbar': tbar }

    mesons = {
        # light unflavoured
        'pip'     : state('pip'  , '$\pi^+$'               , [u, dbar] ),
        'pim'     : state('pim'  , '$\pi^-$'               , [ubar, d] ),
        'piz'     : state('piz'  , '$\pi^0$'               , [u, ubar] ),
        'eta'     : state('eta'  , '$\eta$'                , [d, dbar] ),
        'etapr'   : state('etapr', '$\eta\'$'              , [s, sbar] ),
        'rhop'    : state('rhop' , r'$\rho^{\!+}$'         , [u, dbar] ),
        'rhom'    : state('rhom' , r'$\rho^{\!-}$'         , [ubar, d] ),
        'rhoz'    : state('rhoz' , r'$\rho^0$'             , [u, ubar] ),
        # strange
        'Kp'      : state('Kp'   , '$K^{\!+}$'             , [u, sbar] ),
        'Km'      : state('Km'   , '$K^{\!-}$'             , [ubar, s] ),
        'Kz'      : state('Kz'   , '$K^0$'                 , [d, sbar] ),
        'Kzb'     : state('Kzb'  , r'$\overline{K}^0$'     , [dbar, s] ),
        'omega'   : state('omega', '$\omega$'              , [d, dbar] ),
        'phi'     : state('phi'  , '$\phi\'$'              , [s, sbar] ),
        'Kstp'    : state('Kstp' , '$K^{*\!+}$'            , [u, sbar] ),
        'Kstm'    : state('Kstm' , '$K^{*\!-}$'            , [ubar, s] ),
        'Kstz'    : state('Kstz' , '$K^{*0}$'              , [d, sbar] ),
        'Kstzb'   : state('Kstzb', r'$\overline{K}^{*0}$'  , [dbar, s] ),
        # charm
        'etac'    : state('etac' , '$\eta_c$'              , [c, cbar] ),
        'Dm'      : state('Dm'   , '$D^{\!-}$'             , [cbar, d] ),
        'Dp'      : state('Dp'   , '$D^{\!+}$'             , [c, dbar] ),
        'Dz'      : state('Dz'   , '$D^{0}$'               , [c, ubar] ),
        'Dzb'     : state('Dzb'  , '$\overline{D}^{0}$'    , [cbar, u] ),
        'Dsm'     : state('Dsm'  , '$D_s^{\!-}$'           , [cbar, s] ),
        'Dsp'     : state('Dsp'  , '$D_s^{\!+}$'           , [c, sbar] ),
        'jpsi'    : state('jpsi' , '$J/\psi$'              , [c, cbar] ),
        'Dstm'    : state('Dstm' , '$D^{*\!-}$'            , [cbar, d] ),
        'Dstp'    : state('Dstp' , '$D^{*\!+}$'            , [c, dbar] ),
        'Dstz'    : state('Dstz' , '$D^{*0}$'              , [c, ubar] ),
        'Dstzb'   : state('Dstzb', '$\overline{D}^{*0}$'   , [cbar, u] ),
        'Dsstm'   : state('Dsstm', '$D_s^{*\!-}$'          , [cbar, s] ),
        'Dsstp'   : state('Dsstp', '$D_s^{*\!+}$'          , [c, sbar] ),
    }

    baryons = {
        # unflavoured (spin 1/2)
        'p'        : state('p'         , '$p$'                    , [u,u,d] ),
        'n'        : state('n'         , '$n$'                    , [u,d,d] ),
        # unflavoured (spin 3/2)
        'delm'     : state('delm'      , '$\Delta^-$'             , [d,d,d] ),
        'delz'     : state('delz'      , '$\Delta^0$'             , [u,d,d] ),
        'delp'     : state('delp'      , '$\Delta^+$'             , [u,u,d] ),
        'delpp'    : state('delpp'     , '$\Delta^{+\!\!+}$'      , [u,u,u] ),
        # strange (spin 1/2)
        'sigp'     : state('sigp'      , '$\Sigma^+$'             , [u,u,s] ),
        'sigm'     : state('sigm'      , '$\Sigma^-$'             , [d,d,s] ),
        'sigz'     : state('sigz'      , '$\Sigma^0$'             , [u,d,s] ),
        'lbz'      : state('lbz'       , '$\Lambda^0$'            , [u,d,s] ),
        'xim'      : state('xim'       , '$\Xi^-$'                , [d,s,s] ),
        'xiz'      : state('xiz'       , '$\Xi^0$'                , [u,s,s] ),
        # strange (spin 3/2)
        'sigstp'   : state('sigstp'    , '$\Sigma^{*+}$'          , [u,u,s] ),
        'sigstm'   : state('sigstm'    , '$\Sigma^{*-}$'          , [d,d,s] ),
        'sigstz'   : state('sigstz'    , '$\Sigma^{*0}$'          , [u,d,s] ),
        'xistm'    : state('xistm'     , '$\Xi^{*-}$'             , [d,s,s] ),
        'xistz'    : state('xistz'     , '$\Xi^{*0}$'             , [u,s,s] ),
        'Omega'    : state('Omega'     , '$\Omega^-$'             , [s,s,s] ),
        # charm (spin 1/2)
        'sigcz'    : state('sigcz'     , '$\Sigma_c^0$'           , [d,d,c] ),
        'sigcp'    : state('sigcp'     , '$\Sigma_c^+$'           , [u,d,c] ),
        'sigcpp'   : state('sigcpp'    , '$\Sigma_c^{+\!\!+}$'    , [u,u,c] ),
        'xicz'     : state('xicz'      , '$\Xi_c^0$'              , [d,s,c] ),
        'xicp'     : state('xicp'      , '$\Xi_c^+$'              , [u,s,c] ),
        'lbcp'     : state('lbcp'      , '$\Lambda_c^+$'          , [u,d,c] ),
        'Omegacz'  : state('Omegacz'   , '$\Omega_c^0$'           , [s,s,c] ),
        'Omegaccp' : state('Omegaccp'  , '$\Omega_{cc}^+$'        , [s,c,c] ),
        'xiccp'    : state('xiccp'     , '$\Xi_{cc}^+$'           , [d,c,c] ),
        'xiccpp'   : state('xiccpp'    , '$\Xi_{cc}^{+\!\!+}$'    , [u,c,c] ),
        'Omegaccpp': state('Omegaccpp' , '$\Omega_{cc}^{+\!\!+}$' , [c,c,c] ),
    }

    print_quarks(quarks)
    print_states(mesons , 'Mesons')
    print_states(baryons, 'Baryons')

    # scalar state charm mesons
    ch_mes_s0 = [ 'pip','pim','piz','Kp', 'Km', 'Kz', 'Kzb', 'eta', 'etapr', 'etac', 'Dm','Dp','Dz','Dzb','Dsm','Dsp']
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111,projection='3d')
    plot3d( [ mesons[key] for key in ch_mes_s0 ], ax=ax, mes=0, bar=None)
    fig.tight_layout()
    fig.savefig('pdf/ch_mes_s0.pdf',bbox_inches='tight',pad_inches=-0.4)
    fig.savefig('png/ch_mes_s0.png',bbox_inches='tight',pad_inches=-0.4)

    # vector state charm mesons
    ch_mes_s1 = [ 'rhop','rhom','rhoz','Kstp', 'Kstm', 'Kstz', 'Kstzb', 'omega', 'phi', 'jpsi', 'Dstm','Dstp','Dstz','Dstzb','Dsstm','Dsstp']
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111,projection='3d')
    plot3d( [ mesons[key] for key in ch_mes_s1 ], ax=ax, mes=1, bar=None)
    fig.tight_layout()
    fig.savefig('pdf/ch_mes_s1.pdf',bbox_inches='tight',pad_inches=-0.4)
    fig.savefig('png/ch_mes_s1.png',bbox_inches='tight',pad_inches=-0.4)

    # 1/2 state charm baryons
    ch_bar_s0 = [ 'p', 'n', 'sigp', 'sigm', 'sigz', 'lbz', 'xim', 'xiz', 'sigcz', 'sigcp','sigcpp','xicz','xicp','lbcp','Omegacz', 'xiccp', 'xiccpp', 'Omegaccp' ]
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111,projection='3d')
    plot3d( [ baryons[key] for key in ch_bar_s0 ], ax=ax, mes=None, bar=0)
    fig.tight_layout()
    fig.savefig('pdf/ch_bar_s0.pdf',bbox_inches='tight',pad_inches=-0.4)
    fig.savefig('png/ch_bar_s0.png',bbox_inches='tight',pad_inches=-0.4)

    # 3/2 state charm baryons
    ch_bar_s1 = [ 'delm','delz','delp','delpp','sigp', 'sigm', 'sigz', 'xim', 'xiz', 'Omega','sigcz', 'sigcp','sigcpp','xicz','xicp','Omegacz', 'xiccp', 'xiccpp', 'Omegaccp','Omegaccpp' ]
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111,projection='3d')
    plot3d( [ baryons[key] for key in ch_bar_s1 ], ax=ax, mes=None, bar=1)
    fig.tight_layout()
    fig.savefig('pdf/ch_bar_s1.pdf',bbox_inches='tight',pad_inches=-0.4)
    fig.savefig('png/ch_bar_s1.png',bbox_inches='tight',pad_inches=-0.4)

    # scalar state strange mesons
    sp0_mes = [ 'pip','pim','piz','Kp', 'Km', 'Kz', 'Kzb', 'eta', 'etapr' ]
    fig,ax = plt.subplots(figsize=(6,6))
    plot( [ mesons[key] for key in sp0_mes ] )
    fig.tight_layout()
    fig.savefig('pdf/st_mes_s0.pdf')
    fig.savefig('png/st_mes_s0.png')

    # vector state strange mesons
    sp1_mes = [ 'rhop', 'rhom', 'rhoz', 'omega', 'phi', 'Kstp', 'Kstm', 'Kstz', 'Kstzb' ]
    fig,ax = plt.subplots(figsize=(6,6))
    plot( [ mesons[key] for key in sp1_mes ] )
    fig.tight_layout()
    fig.savefig('pdf/st_mes_s1.pdf')
    fig.savefig('png/st_mes_s1.png')

    # 1/2 state strange baryons
    sc_bar = [ 'p', 'n', 'sigp', 'sigm', 'sigz', 'lbz', 'xim', 'xiz' ]
    fig, ax = plt.subplots(figsize=(6,6))
    plot( [ baryons[key] for key in sc_bar ], hexc=(0,-1) )
    fig.tight_layout()
    fig.savefig('pdf/st_bar_s0.pdf')
    fig.savefig('png/st_bar_s0.png')

    # 3/2 state strange baryons
    sv_bar = [ 'delm', 'delz', 'delp', 'delpp', 'sigstm', 'sigstz', 'sigstp', 'xistm', 'xistz', 'Omega' ]
    fig, ax = plt.subplots(figsize=(6,6))
    plot( [ baryons[key] for key in sv_bar ], hexc=(0,-1), hexr=None, tri=((-1.5,1.5,0,-1.5),(0,0,-3,0)) )
    fig.tight_layout()
    fig.savefig('pdf/st_bar_s1.pdf')
    fig.savefig('png/st_bar_s1.png')

    # show plots if -i passed
    if len(sys.argv)>1:
        if sys.argv[1] == '-i':
            plt.show()
