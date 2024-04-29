import time
from pfac.fac import *

file_name = 'o8_he'

nele = 1

if nele==1:
    config_labels_lower = ['n0']
    config_labels_upper = ['n1','n2','n3','n4','n5','n6','n7','n8','n9']

elif nele==2:
    config_labels_lower = ['n1']
    config_labels_upper = ['n11','n12','n13','n14','n15','n16','n17','n18','n19']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save_dir = f'/home/tim/research/cuebit_col/fac_output/o8_he/'

t = time.time()
t1 = t
file = save_dir+file_name

def e_configs(config_labels):
    econf = []
    for i in range(len(config_labels)):
        econf.append('')
        for n in range(1,10):
            is_e = config_labels[i].count(str(n))
            if is_e != 0:
                econf[i] = econf[i] + str(n) + '*' + str(is_e) + ' '
        econf[i] = econf[i].strip()
    return econf

def lap(label):
    global t1
    print(f'{label}: {time.time()-t1:.2f} s')
    t1 = time.time()

InitializeMPI(22)
SetAtom('O')

for label, econf in zip(config_labels_lower, e_configs(config_labels_lower)):
    Config(label, econf)

for label, econf in zip(config_labels_upper, e_configs(config_labels_upper)):
    Config(label, econf)

lap('Init')

#ConfigEnergy(0)
OptimizeRadial(config_labels_upper[0])
#ConfigEnergy(1)
Structure(file+'.b.en', config_labels_upper)
MemENTable(file+'.b.en')
PrintTable(file+'.b.en',file+'.en',1)

lap('Upper Struct')

#ConfigEnergy(0)
#OptimizeRadial(config_labels_lower[0])
#ConfigEnergy(1)
Structure(file+'.b.en', config_labels_lower)
MemENTable(file+'.b.en')
PrintTable(file+'.b.en',file+'.en',1)

lap('Lower Struct')

TransitionTable(file+'.b.tr', config_labels_upper, config_labels_upper)
PrintTable(file+'.b.tr',file+'.tr', 1)

lap('Trans')

SetCXTarget('He')
SetCXEGrid(100,1000,100000)
CXTable(file+'.b.cx', config_labels_upper, config_labels_lower)
PrintTable(file+'.b.cx',file+'.cx', 1)

lap('CX')

FinalizeMPI()
print(f'FAC Total: {time.time()-t:.3f}')