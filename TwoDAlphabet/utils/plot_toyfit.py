'''Plot a toy data set against model from which it was generated.

Highly repetetive of the 2D Alphabet code base and should be re-written after refactor.
'''

import header, os, pickle
import ROOT
from ROOT import *
from optparse import OptionParser
from array import array

gStyle.SetOptStat(0)
gROOT.SetBatch(kTRUE)

parser = OptionParser()
parser.add_option("-p", "--projPath", dest="projPath",
                help="Home of the project - has the cards, fit results, etc")
parser.add_option('-r', '--regions', dest='regions',
                help='Comma separated list of names of regions fit from configs (name in each json config)')
parser.add_option('-n', '--name', dest='name',
                help='Naming of output files')
parser.add_option('-w', '--workspacefile', dest='workspacefile', default='',
                help='Naming of workspace file')
parser.add_option('--toyfile', dest='toyfile', default='',
                help='Naming of toy file')
parser.add_option('-s', '--seed', dest='seed',
                help='Seed')
parser.add_option('-t', '--toynum', dest='toynum',
                help='Toy number')
parser.add_option("--skipSamples", action="store_true", 
                default =   False,
                dest    =   "skipSamples",
                help    =   "Skip PostFit2DShapesFromWorkspace step")
(options, args) = parser.parse_args()

class ToyFit():
    def __init__ (self,options):
        self.projPath = options.projPath+'/'
        self.regions = options.regions
        self.name = options.name
        self.toyDir = 'toy'+options.toynum+'_'+self.name+'_plots/'
        self.toynum = options.toynum
        if not os.path.isdir(self.projPath+'/'+self.toyDir): os.mkdir(self.projPath+'/'+self.toyDir)

        
        if os.path.isfile(self.projPath+'/fitDiagnostics'+self.name+'.root'): self.fd_file = TFile.Open(self.projPath+'/fitDiagnostics'+self.name+'.root')
        else: self.fd_file = False
        if options.workspacefile == '': self.ws_file = TFile.Open(self.projPath+'/higgsCombine'+self.name+'.GoodnessOfFit.mH120.'+options.seed+'.root')
        else: self.ws_file = TFile.Open(self.projPath+'/'+options.workspacefile)
        if options.toyfile == '': self.toy_file = self.ws_file
        else: self.toy_file = TFile.Open(self.projPath+'/'+options.toyfile)
        # self.ws_file.ls()
        self.ws = self.ws_file.Get('w')
        # self.ws.Print()
        # self.toy_w = self.ws.Clone('toy_w')
        self.toy = self.toy_file.Get('toys/toy_'+options.toynum)
        getattr(self.ws,'import')(self.toy)
        self.ws.writeToFile(self.projPath+'/'+self.toyDir+'/toy_workspace.root',True)

        if not options.skipSamples:
            with header.cd(self.projPath):
                header.executeCmd('PostFit2DShapesFromWorkspace -w '+self.toyDir+'/toy_workspace.root --dataset model_sData -o '+self.toyDir+'postfitshapes_s.root -f fitDiagnostics'+self.name+'.root:fit_s --samples 50 --postfit --sampling --print')
        with header.cd(self.projPath):
            # diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics'+self.name+'.root --abs -g toy'+options.toynum+'_'+self.name+'_plots/nuisance_pulls.root'
            # header.executeCmd(diffnuis_cmd)
            self.plotNuisances()


        self.post_file = TFile.Open(self.projPath+self.toyDir+'postfitshapes_s.root')
        for r in self.regions.split(','):
            if not os.path.isdir(self.projPath+'/'+self.toyDir+r): os.mkdir(self.projPath+'/'+self.toyDir+r)
            if os.path.isfile(self.projPath+'/saveOut.p'): pickleFile = pickle.load(open(self.projPath+'/saveOut.p','rb'))
            else: pickleFile = pickle.load(open(self.projPath+r+'/saveOut.p','rb'))
            self.plotFitResults(r,pickleFile)

    def plotNuisances(self):
        nuis_list = []
        tree_fit_sb = self.fd_file.Get('tree_fit_sb')
        tree_fit_sb.GetEntry(int(self.toynum)-1)
        branches = tree_fit_sb.GetListOfBranches()
        for ibranch in range(branches.GetEntries()):
            branch_name = branches.At(ibranch).GetName()
            if '_In' in branch_name:
                print ('Adding nuisance: '+branch_name.replace('_In',''))
                nuis_list.append(branch_name.replace('_In',''))

        out_plot = TH1F('nuisance_pulls','nuisance_pulls',len(nuis_list),0,len(nuis_list))

        for i,n in enumerate(nuis_list):
            out_plot.SetBinContent(i+1,getattr(tree_fit_sb,n))
            out_plot.SetBinError(i+1,1)
            out_plot.GetXaxis().SetBinLabel(i+1,n)

        c_out = TCanvas('nuisance_pulls','Nuisance Parameters',1000,600)
        out_plot.SetMinimum(-4)
        out_plot.SetMaximum(4)
        out_plot.Draw('pe')
        c_out.Print('toy'+options.toynum+'_'+self.name+'_plots/nuisance_pulls.pdf','pdf')

    def plotFitResults(self,region,pickle_file): # fittag means 'b' or 's'
        allVars = []
        fittag = 's'
        plotOutDir = self.projPath+self.toyDir+region
        #####################
        #   Get everything  #
        #####################

        # File with histograms and RooFitResult parameters
        x_low = pickle_file['newXbins']['LOW'][0]
        x_high = pickle_file['newXbins']['HIGH'][-1]
        y_low = pickle_file['newYbins'][0]
        y_high = pickle_file['newYbins'][-1]
        y_nbins = len(pickle_file['newYbins'])-1

        fullXbins = pickle_file['full_x_bins']
        newYbins = pickle_file['newYbins']
        sigStart = pickle_file['sigStart']
        sigEnd = pickle_file['sigEnd']
        xVarTitle = pickle_file['xVarTitle']
        yVarTitle = pickle_file['yVarTitle']

        name = pickle_file['name']

        if os.path.isfile(self.projPath+'/saveOut.p'): inputConfig = header.openJSON(self.projPath+'/runConfig.json')
        else: inputConfig = header.openJSON(self.projPath+region+'/runConfig.json')
        blindedPlots = False#inputConfig['OPTIONS']['blindedPlots']

        # Define low, middle, high projection regions for y (x regions defined already via signal region bounds)
        y_turnon_endBin = self.post_file.Get('pass_LOW_'+name+'_prefit/data_obs').ProjectionY().GetMaximumBin()
        y_tail_beginningBin = int((y_nbins - y_turnon_endBin)/2.0 + y_turnon_endBin)
        print ('Finding start and end bin indexes of signal range. Looking for '+str(sigStart)+', '+str(sigEnd))
        for ix,xwall in enumerate(fullXbins):
            if xwall == sigStart:
                print ('Assigning start bin as '+str(ix+1))
                x_sigstart_bin = ix+1
            if xwall == sigEnd:
                print ('Assigning end bin as '+str(ix))
                x_sigend_bin = ix

        if y_turnon_endBin > y_nbins/2.0:  # in case this isn't a distribution with a turn-on
            y_turnon_endBin = int(round(y_nbins/3.0))
            y_tail_beginningBin = 2*y_turnon_endBin
        y_turnon_endVal = str(int(self.post_file.Get('pass_LOW_'+name+'_prefit/data_obs').GetYaxis().GetBinUpEdge(y_turnon_endBin)))
        y_tail_beginningVal = str(int(self.post_file.Get('pass_LOW_'+name+'_prefit/data_obs').GetYaxis().GetBinLowEdge(y_tail_beginningBin)))
     

        # Final fit signal strength
        if self.fd_file:
            if fittag == 's':
                tree_fit_sb = self.fd_file.Get('tree_fit_sb')
                tree_fit_sb.GetEntry(int(self.toynum)-1)
                signal_strength = tree_fit_sb.r
            else:
                tree_fit_b = self.fd_file.Get('tree_fit_b')
                tree_fit_b.GetEntry(int(self.toynum)-1)
                signal_strength = tree_fit_b.r
        else: 
            signal_strength = 0

        #####################
        #    Data vs Bkg    #
        #####################

        hist_dict = {}

        # Organize and make any projections or 2D distributions
        for process in [process for process in inputConfig['PROCESS'] if process != 'HELP']+['qcd']:
            hist_dict[process] = {}
            for cat in ['fail','pass']:
                hist_dict[process][cat] = {'LOW':{},'SIG':{},'HIGH':{}}
                x_slice_list_pre = []
                x_slice_list_post = []
                # Grab everything and put clones in a dictionary
                for c in ['LOW','SIG','HIGH']:
                    file_dir = cat+'_'+c+'_'+name
                    hist_dict[process][cat]['prefit_'+c] = self.post_file.Get(file_dir+'_prefit/'+process).Clone()
                    hist_dict[process][cat]['postfit_'+c] = self.post_file.Get(file_dir+'_postfit/'+process).Clone()
                    x_slice_list_pre.append(hist_dict[process][cat]['prefit_'+c])    # lists for 2D making
                    x_slice_list_post.append(hist_dict[process][cat]['postfit_'+c])

                # First rebuild the 2D distributions
                if blindedPlots and process == 'data_obs':
                    hist_dict[process][cat]['prefit_2D'] = header.stitchHistsInX(process+'_'+cat+'_prefit2D',fullXbins,newYbins,x_slice_list_pre,blinded=[1])
                    hist_dict[process][cat]['postfit_2D'] = header.stitchHistsInX(process+'_'+cat+'_postfit2D',fullXbins,newYbins,x_slice_list_post,blinded=[1])

                else:
                    hist_dict[process][cat]['prefit_2D'] = header.stitchHistsInX(process+'_'+cat+'_prefit2D',fullXbins,newYbins,x_slice_list_pre,blinded=[])
                    hist_dict[process][cat]['postfit_2D'] = header.stitchHistsInX(process+'_'+cat+'_postfit2D',fullXbins,newYbins,x_slice_list_post,blinded=[])

                hist_dict[process][cat]['prefit_2D'].SetMinimum(0)
                hist_dict[process][cat]['postfit_2D'].SetMinimum(0)
                hist_dict[process][cat]['prefit_2D'].SetTitle(process + ', ' + cat +', '+name+ ', pre-fit')
                hist_dict[process][cat]['postfit_2D'].SetTitle(process + ', ' + cat +', '+name+ ', post-fit')

                # Now projections
                base_proj_name_pre = process+'_'+cat+'_'+name+'_pre_'
                base_proj_name_post = process+'_'+cat+'_'+name+'_post_'

                hist_dict[process][cat]['prefit_projx1'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+str(y_low)+'-'+y_turnon_endVal,              1,                      y_turnon_endBin, 'e')
                hist_dict[process][cat]['prefit_projx2'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,     y_turnon_endBin+1,      y_tail_beginningBin,'e')
                hist_dict[process][cat]['prefit_projx3'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_tail_beginningVal+'-'+str(y_high),         y_tail_beginningBin+1,  y_nbins,'e')

                hist_dict[process][cat]['prefit_projy1'] = hist_dict[process][cat]['prefit_LOW'].ProjectionY(base_proj_name_pre+'projy_'+str(x_low)+'-'+str(sigStart),          1,                      hist_dict[process][cat]['prefit_LOW'].GetNbinsX(),'e')
                if blindedPlots:
                    hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(sigStart)+'-'+str(sigEnd), x_sigstart_bin,           x_sigend_bin,'e')
                else:
                    hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_SIG'].ProjectionY(base_proj_name_pre+'projy_'+str(sigStart)+'-'+str(sigEnd),    1,                      hist_dict[process][cat]['prefit_SIG'].GetNbinsX(),'e')
                hist_dict[process][cat]['prefit_projy3'] = hist_dict[process][cat]['prefit_HIGH'].ProjectionY(base_proj_name_pre+'projy_'+str(sigEnd)+'-'+str(x_high),          1,                      hist_dict[process][cat]['prefit_HIGH'].GetNbinsX(),'e')

                hist_dict[process][cat]['postfit_projx1'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+str(y_low)+'-'+y_turnon_endVal,           1,                      y_turnon_endBin, 'e')
                hist_dict[process][cat]['postfit_projx2'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,  y_turnon_endBin+1,      y_tail_beginningBin,'e')
                hist_dict[process][cat]['postfit_projx3'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_tail_beginningVal+'-'+str(y_high),      y_tail_beginningBin+1,  y_nbins,'e')

                hist_dict[process][cat]['postfit_projy1'] = hist_dict[process][cat]['postfit_LOW'].ProjectionY(base_proj_name_post+'projy_'+str(x_low)+'-'+str(sigStart),       1,                      hist_dict[process][cat]['postfit_LOW'].GetNbinsX(),'e')
                if blindedPlots:
                    hist_dict[process][cat]['postfit_projy2'] = hist_dict[process][cat]['postfit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(sigStart)+'-'+str(sigEnd), x_sigstart_bin,           x_sigend_bin,'e')
                else:
                    hist_dict[process][cat]['postfit_projy2'] = hist_dict[process][cat]['postfit_SIG'].ProjectionY(base_proj_name_post+'projy_'+str(sigStart)+'-'+str(sigEnd), 1,                      hist_dict[process][cat]['postfit_SIG'].GetNbinsX(),'e')
                hist_dict[process][cat]['postfit_projy3'] = hist_dict[process][cat]['postfit_HIGH'].ProjectionY(base_proj_name_post+'projy_'+str(sigEnd)+'-'+str(x_high),       1,                      hist_dict[process][cat]['postfit_HIGH'].GetNbinsX(),'e')

                x_edges = [x_low,sigStart,sigEnd,x_high]
                y_edges = [y_low,y_turnon_endVal,y_tail_beginningVal,y_high]

                for z in ['x','y']:
                    for i in range(1,4):
                        hist_dict[process][cat]['postfit_proj'+z+str(i)].SetMinimum(0)
                        if z == 'x':
                            hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+', '+name+ ', ' +str(y_edges[i-1]) +'-'+ str(y_edges[i]))
                        elif z == 'y':
                            hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+', '+name+ ', ' +str(x_edges[i-1]) +'-'+ str(x_edges[i]))

        #post_file.Close()

        # Add together processes that we want to see as one
        # if self.plotTogether != False:
        #     hist_dict = self.plotProcessesTogether(hist_dict)
            
        process_list = hist_dict.keys()

        # Create lists for the 2D projections (ones you want to see together)
        for process in hist_dict.keys():    # Canvas
            isSignal = (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] == 0)
            twoDList = []         
            for cat in ['fail','pass']:
                for fit in ['prefit', 'postfit']:
                    if isSignal and fittag == 's' and fit == 'postfit':
                        # hist_dict[process][cat][fit+'_2D'].Scale(signal_strength)
                        twoDList.append(hist_dict[process][cat][fit+'_2D'])
                    else:
                        twoDList.append(hist_dict[process][cat][fit+'_2D'])

            if isSignal and fittag != 's':
                continue
            else:
                header.makeCan(process+'_fit'+fittag+'_2D',plotOutDir,twoDList,xtitle=xVarTitle,ytitle=yVarTitle)

        # Invert the last two items (unique to b*) - customize as needed
        process_list[-1],process_list[-2] = process_list[-2],process_list[-1]

        # Get the colors
        colors = []
        for process in process_list:
            if process != 'data_obs':
                if process not in inputConfig['PROCESS'].keys():
                    continue
                if (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] != 0):
                    if (process in inputConfig['PROCESS'].keys()) and ('COLOR' in inputConfig['PROCESS'][process].keys()):
                        colors.append(inputConfig['PROCESS'][process]["COLOR"])
                    else:
                        colors.append(None)

        # Put QCD on bottom of stack since it's smooth
        colors = [kYellow]+colors

        # Create lists for makeCan of the projections
        for plotType in ['postfit_projx','postfit_projy']:   # Canvases
            bkgList = []
            dataList = []
            signal_list = []
            
            for cat in ['fail','pass']: # Row 
                for regionNum in range(1,4):    # Column
                    bkg_process_list = []
                    for process in process_list:
                        if process != 'data_obs':
                            if (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] != 0):
                                bkg_process_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                            elif (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] == 0):
                                # if fittag == 'fit_s':
                                # hist_dict[process][cat][plotType+str(regionNum)].Scale(signal_strength)
                                signal_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                                # else:
                                    # signal_list.append(hist_dict[process][cat][plotType+str(regionNum)]) #.replace('post','pre')

                        else:
                            dataList.append(hist_dict[process][cat][plotType+str(regionNum)])

                    
                    ## Need to do this ##
                    # if 'x' in plotType:
                    #     hist_dict['qcd'][cat][plotType+str(regionNum)].GetXaxis().Set(len(self.fullXbins)-1,array.array('d',self.fullXbins))
                    # elif 'y' in plotType:
                    #     hist_dict['qcd'][cat][plotType+str(regionNum)].GetXaxis().Set(len(self.newYbins)-1,array.array('d',self.newYbins))
                    ## otherwise a bunch of dumb warnings are thrown about bin limits ##

                    # Put QCD on bottom of stack since it's smooth
                    bkg_process_list = [hist_dict['qcd'][cat][plotType+str(regionNum)]]+bkg_process_list


                    bkgList.append(bkg_process_list)

            # An attempt to debug the warning about adding histograms with different bin limits - seems there's no reason the warning should be coming up
            # for i,data in enumerate(dataList):
            #     print data.GetName()
            #     for b in bkgList[i]:
            #         for i in range(data.GetXaxis().GetXbins().fN):
            #             if not TMath.AreEqualRel(data.GetXaxis().GetXbins().GetAt(i),b.GetXaxis().GetXbins().GetAt(i),0.00000000001):
            #                 print b.GetName() + ' does not agree. Data: '+str(data.GetXaxis().GetXbins().GetAt(i)) +' vs Bkg: '+str(b.GetXaxis().GetXbins().GetAt(i))


            if 'x' in plotType:
                header.makeCan(plotType+'_fit'+fittag,plotOutDir,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=xVarTitle)
                header.makeCan(plotType+'_fit'+fittag+'_log',plotOutDir,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=xVarTitle,logy=True)
            elif 'y' in plotType:
                header.makeCan(plotType+'_fit'+fittag,plotOutDir,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=yVarTitle)
                header.makeCan(plotType+'_fit'+fittag+'_log',plotOutDir,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=yVarTitle,logy=True)

        # Make comparisons for each background process of pre and post fit projections
        for plotType in ['projx','projy']:
            for process in process_list:
                if process != 'data_obs':
                    pre_list = []
                    post_list = []
                    for cat in ['fail','pass']: # Row 
                        for regionNum in range(1,4):    # Column
                            pre_list.append([hist_dict[process][cat]['prefit_'+plotType+str(regionNum)]])  # in terms of makeCan these are "bkg hists"
                            post_list.append(hist_dict[process][cat]['postfit_'+plotType+str(regionNum)])   # and these are "data hists"
                            if process != 'qcd':
                                if 'COLOR' in inputConfig['PROCESS'][process].keys():
                                    prepostcolors = [inputConfig['PROCESS'][process]['COLOR']]
                                else:
                                    prepostcolors = [0]
                            else:
                                prepostcolors = [kYellow]

                    if 'x' in plotType: header.makeCan(process+'_'+plotType+'_fit'+fittag,plotOutDir,post_list,bkglist=pre_list,colors=prepostcolors,xtitle=xVarTitle,datastyle='histe')
                    if 'y' in plotType: header.makeCan(process+'_'+plotType+'_fit'+fittag,plotOutDir,post_list,bkglist=pre_list,colors=prepostcolors,xtitle=yVarTitle,datastyle='histe')

print (options)
tf = ToyFit(options)
