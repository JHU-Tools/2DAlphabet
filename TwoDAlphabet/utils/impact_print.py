'''Print out nuisance parameter impacts as a table instead of a plot.
Averages impacts across years if possible. System argument is the output 
json from the combineTool.py Impacts method.
'''

import sys
sys.path.append('../')
import header

f = header.openJSON(sys.argv[1],False)

param_names_all = {}
param_names = {}

for pdict in f['params']:
    pname_noyear = pdict['name'].replace('16','').replace('17','').replace('18','')
    param_names_all[pdict['name']] = {"up":
                                        (pdict['r'][0]-pdict['r'][1])/pdict['r'][1],
                                      "down":
                                        (pdict['r'][2]-pdict['r'][1])/pdict['r'][1]}
    if pname_noyear not in param_names.keys():
        param_names[pname_noyear] = {'up':0,'down':0}

for p in param_names.keys():
    n = 0
    for k in param_names_all.keys():
        if p in k:
            print ('\t +%s(%.2f,%.2f)'%(k, param_names_all[k]['up']*100, param_names_all[k]['down']*100) )
            n+=1
            param_names[p]['up'] += param_names_all[k]['up']
            param_names[p]['down'] += param_names_all[k]['down']
    if n > 1: print ('Found %s %s'%(n,p))
    param_names[p]['up'] /= n
    param_names[p]['down'] /= n

for p in sorted(param_names.keys(), key=str.lower):
    print ('%s: %.1f/%.1f'%(p,param_names[p]['up']*100,param_names[p]['down']*100))

print ('========ALL========')
for p in sorted(param_names_all.keys(), key=str.lower):
    print ('%s: %.1f/%.1f'%(p,param_names_all[p]['up']*100,param_names_all[p]['down']*100))
