from datatools import *
from plot_benchcore import *
fin=open("../Data/EVTS.dat",'r')


dico={}
test=10000000
for l in fin:
    line=l.split()
    if not line:
        continue
    typ=line[0]
    tup=tuple(eval(line[1]))
    dico.setdefault(typ,{})
    dico[typ].setdefault(tup,0)
    dico[typ][tup]+=1
    if test<0:
        break
    test-=1
fin.close()

#Number of occupied schools
no=[]
for typ,dic in dico.iteritems():
    for t in dic:
        no.append((t[0]*t[2],t[1]*t[2],t[4],t[5],t[6]))
no=array(no).T
#scatter(no[1]*(1-np.exp(-no[0]/no[1])),no[2],hold=1)
ns=no[2]*(1-np.exp(-no[0]/no[2]))
plot(*make_density((ns,no[3]),bins=20)[:2],xlabel='n_o (theo)',ylabel='n_o (exp)')
scatter((no[1]+no[0])[::10],no[-1][::10])
#plot(*make_density(( meano[meanox==False],no[-1][meanox==False]),bins=20)[:2],xlabel='<o> (theo)',ylabel='<o> (exp)')


for typ,dic in dico.iteritems():
    sb=array([(t[0]*t[2],t[1]*t[2],dic[t]) for t in dic][::40],dtype="float")
    s={t[0]*t[2]:dic[t] for t in dic}
    b={t[1]*t[2]:dic[t] for t in dic}
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel("finder")
    plt.ylabel("bound")
    ax.set_zlabel(typ)
    ax.scatter(sb[:,0],sb[:,1],sb[:,2])
    plt.show()
    x,y=make_density(list(s.iteritems()))[:2]
    fout=open("parsed_{}.dat".format(typ),'w')
    for i in range(len(x)):
        fout.write("{} {}\n".format(x,y))
    fout.close()
    plt.title(typ)
    plot(x,y)
