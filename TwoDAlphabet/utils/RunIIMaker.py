import ROOT, sys, pprint
from ROOT import TFile
pp = pprint.PrettyPrinter(indent = 2)

def RunIIMaker(projDirs):
    for projDir in projDirs:
        for fittype in ['b','s']:
            plot_file = TFile.Open(projDir+'/postfitshapes_%s.root'%fittype)
            file_dirs = plot_file.GetListOfKeys()
            out_file = TFile.Open(projDir+'/postfitshapes_RunII_%s.root'%fittype,'RECREATE')

            # Get tags from directory structure of TFile
            structure = {}
            for prepost in ['_prefit','_postfit']:
                for k in file_dirs:
                    if prepost in k.GetName():
                        # dir names formated as something like fail_HIGH_<tag we want>_<prepost>
                        category = '_'.join(k.GetName().split('_')[0:2])
                        this_tag = ''.join(k.GetName().split('_')[2:-1])

                        if category not in structure.keys(): 
                            structure[category] = {}

                        if ('2016' in this_tag or '2017' in this_tag or '2018' in this_tag):
                            this_tag = this_tag.replace('2016','16').replace('2017','17').replace('2018','18')

                        if '16' in this_tag: year = '16'
                        elif '17' in this_tag: year = '17'
                        elif '18' in this_tag: year = '18'
                        else: year = ''

                        if year != '':
                            this_tag = this_tag.replace(year,'') # drop any possible year

                        if this_tag not in structure[category].keys(): 
                            structure[category][this_tag] = {}

                        for proc in plot_file.Get(k.GetName()).GetListOfKeys():
                            if proc.GetName() not in structure[category][this_tag].keys():
                                structure[category][this_tag][proc.GetName()] = {}

                            structure[category][this_tag][proc.GetName()][year] = plot_file.Get(k.GetName()+'/'+proc.GetName())
                            
                for category in structure.keys():
                    for tag in structure[category].keys():
                        for proc in structure[category][tag].keys():
                            years = structure[category][tag][proc].keys()
                            base = structure[category][tag][proc][years[0]].Clone()
                            for year in years[1:]:
                                base.Add(structure[category][tag][proc][year])

                            base.SetTitle(base.GetTitle().replace(years[0],'RunII'))
                            new_dir_name = '%s_%s%s'%(category,tag+'RunII',prepost)
                            if new_dir_name not in [k.GetName() for k in out_file.GetListOfKeys()]:
                                out_file.mkdir(new_dir_name)
                            out_file.cd(new_dir_name)
                            # base.SetDirectory(ROOT.gDirectory)
                            base.Write()

            out_file.Close()

        # return structure[category].keys()

if __name__ == '__main__':
    RunIIMaker(sys.argv[1:])
