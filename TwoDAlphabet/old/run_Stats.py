import sys, os, time, decimal, pickle
import header
import random
import ROOT
from ROOT import *

gStyle.SetOptStat(0)
gStyle.SetOptFit(1)
gStyle.SetOptTitle(0)
gROOT.SetBatch()

def getMasks(filename):
    # Determine if there are channel masks from the fit in the projDir
    masked_regions = []
    f = TFile.Open(filename)
    w = f.Get('w')
    allVars = RooArgList(w.allVars())

    # Loop over all vars in original workspace and add any masked region names to list
    for i in range(allVars.getSize()):
        if 'mask_pass_SIG_' in allVars[i].GetName():
            if allVars[i].getValV() == 1:
                masked_regions.append(allVars[i].GetName())

    f.Close()

    return masked_regions
   

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--projDir", metavar='<dir>', dest="projDir",default='',
                help="Home of the project - has the cards, fit results, etc")
parser.add_option("-a", "--altDir", metavar='<dir>', dest="altDir",default='',
                help="Home of the alternative model that you'd like to compare against the one in projDir")
parser.add_option("-t", "--toys", metavar='<N>', dest="toys",default='100',
                help="Number of toys to generate - fed to combine")
parser.add_option("--toyJobs", metavar='<N>', dest="toyjobs",default='1',
                help="Number of jobs to split toy fitting into - fed to condor job building")
parser.add_option("--seed", metavar='<N>', dest="seed",default='',
                help="Seed - fed to combine")
parser.add_option('--rMin', metavar='<rMin>', type='string', action='store',
                default =   '-5',
                dest    =   'rMin',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option('--rMax', metavar='<rMax>', type='string', action='store',
                default =   '5',
                dest    =   'rMax',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option('-w', metavar='<name>', type='string', action='store',
                default =   '',
                dest    =   'workspace',
                help    =   'Model workspace (from text2workspace)')
# parser.add_option("-b", "--blind", 
#                 action="store_false", dest="blind", default=True,
#                 help="Blinds the signal region in the fit")
parser.add_option("--plotOnly", 
                action="store_true", dest="plotOnly", default=False,
                help="Only plot")
parser.add_option("--dryrun", 
                action="store_true", dest="dryrun", default=False,
                help="Dry run the combine commands to console")
parser.add_option("--gof", 
                action="store_true", dest="gof", default=False,
                help="Perform goodness of fit test")
parser.add_option("--signalInjection", 
                action="store", metavar='<r>', dest="signalInjection", default='',
                help="Perform signal injection test")
parser.add_option("--biasStudy", 
                action="store", dest="biasStudy", default='',
                help="Perform bias study")
parser.add_option("--ftest", 
                action="store", metavar='<option>', dest="ftest", default=False,
                help="Perform F test. Options are 'generate', 'fitAlt', 'fitMain', 'post', 'pvalue'")
# parser.add_option("--diagnosticsWithToys", 
#                 action="store_true", dest="diagnosticsWithToys", default=False,
#                 help="Perform diagnostics with toys")
parser.add_option("--post", 
                action="store_true", dest="post", default=False,
                help="Run in conjunction with diagnosticsWithToys or signalInjection and condor jobs to process output files")
parser.add_option("--freezeFail", 
                action="store_true", dest="freezeFail", default=False,
                help="Run in conjunction with diagnosticsWithToys or signalInjection and condor jobs to process output files")
parser.add_option("--condor", 
                action="store_true", dest="condor", default=False,
                help="Submit a condor job (only for signal injection and diagnosticsWithToys currently")
parser.add_option("--skipSnapshot", 
                action="store_true", dest="skipSnapshot", default=False,
                help="Use if you've already made a snapshot that you trust for this project directory")
(options, args) = parser.parse_args()

projDir = options.projDir # home of the workspace - has the cards, fit results, etc
tag = projDir.split('/')[0]

if projDir.split('/')[-1] == '': card_tag = projDir.split('/')[-2]
else: card_tag = projDir.split('/')[-1]

if tag == '':
    print 'ERROR in project directory name (where your workspace and data card lives). Did you accidentally provide a leading slash? (ie /projDir/) Quitting...'
    quit()

if not os.path.isdir(projDir): 
    print projDir +' is not a directory. Quitting...'
    quit()

if (options.biasStudy != '' or options.ftest == 'fitAlt' or options.ftest == 'post' or options.ftest == 'pvalue'):
    altDir = options.altDir # home of the alternate workspace - has the cards, fit results, etc
    alttag = altDir.split('/')[0]
    if altDir.split('/')[-1] != '': altcard_tag = altDir.split('/')[-1]
    else: altcard_tag = altDir.split('/')[-2]

    if altDir[-1] != '/': altDir_depth = '../'*(altDir.count('/')+1)
    else: altDir_depth = '../'*(altDir.count('/'))

    if not os.path.isdir(altDir): 
        print altDir +' is not a directory. Quitting...'
        quit()

# Start with tests that only require the projDir and don't do a comparison
with header.cd(projDir):
    # Determine if there are channel masks from the fit in the projDir
    masked_regions = getMasks('higgsCombineTest.FitDiagnostics.mH120.root')
    mask_string = '' if len(masked_regions) == 0 else ' --setParameters '
    for mask in masked_regions: mask_string+=mask+'=1,'
    mask_string = mask_string[:-1]

    # Set the string to specify the blinding
    blind_string = "--setParametersForFit "
    for m in masked_regions:
        blind_string+=m+'=1,'
    blind_string = blind_string[:-1]+' '
    blind_string = blind_string+blind_string.replace('setParametersForFit','setParametersForEval')#.replace('=1','=0')

    if len(masked_regions) > 0: 
    	freeze_r_string = " --fixedSignalStrength 0"
    else: freeze_r_string = ''

    if options.freezeFail:
        freeze_string = ' --freezeParameters "var{Fail_.*}"'
    else: freeze_string = ''

    # Make a prefit workspace from the data card
    if options.workspace == '' and not options.post:
        workspace_name = 'stats_workspace.root'
        t2w_cmd = 'text2workspace.py --channel-masks -b card_'+card_tag+'.txt -o '+workspace_name
        # if not (os.path.exists(workspace_name)):# and options.ftest == 'pvalue'):
        header.executeCmd(t2w_cmd,options.dryrun)
    else:
        workspace_name = options.workspace

    # Setup a random seed
    seed = random.randint(100000,999999)
    if options.seed != '': seed = options.seed

    # Setup run naming
    if options.gof:
        expectSignal = '0'
        run_name = 'gof'
        gen_name = run_name
    elif options.signalInjection != '':
        expectSignal = options.signalInjection
        run_name = 'signalInjection'+expectSignal.replace('.','p')
        gen_name = run_name
    elif options.biasStudy != '':
        expectSignal = options.biasStudy
        run_name = 'biasStudy'+expectSignal.replace('.','p')+'_'+tag+'v'+alttag
        gen_name = run_name
    elif options.ftest:
        expectSignal = '0'
        run_name = 'FTest_'+tag
        gen_name = 'FTest'

    # Setup toy naming and values and command to run them
    run_name+='_'+options.toys
    ntoys = int(options.toys)
    toyjobs = int(options.toyjobs) if not options.ftest else 1
    toys_per_job = int(ntoys/toyjobs)
    toy_dict = {
        'ntoys':ntoys,
        'toyjobs':toyjobs,
        'toys_per_job':toys_per_job,
        'seed':seed
    }
    gen_command = 'combine -M GenerateOnly'+mask_string.replace('=1','=0')+' -d initialFitWorkspace.root --snapshotName initialFit --toysFrequentist --bypassFrequentistFit -t '+str(ntoys)+' --saveToys -s '+str(seed)+' --expectSignal '+expectSignal+' -n '+gen_name+freeze_string

    if not options.post and not options.plotOnly and not options.skipSnapshot: header.setSnapshot()

    #######
    # GOF #
    #######
    if options.gof:
        if not options.plotOnly and not options.post:
            commands = []

            ########################################################################
            # First morph the base workspace to post-fit according to MLfit result #
            ########################################################################
            # Check if we can import post-fit result made during MLfit step
            # if not os.path.isfile('fitDiagnosticsTest.root'):
            #     print 'ERROR: '+projDir+'/fitDiagnosticsTest.root does not exist. Please check that run_MLfit.py finished correctly. Quitting...'
            #     quit()

            gof_data_cmd = 'combine -M GoodnessOfFit '+workspace_name+' --algo=saturated '+blind_string+freeze_string+freeze_r_string+' -n gof_data'
            header.executeCmd(gof_data_cmd,options.dryrun)
            gof_toy_cmd = 'combine -M GoodnessOfFit '+workspace_name+' --algo=saturated --toysFrequentist --toysFile higgsCombine'+gen_name+'.GenerateOnly.mH120.'+str(seed)+'.root --saveWorkspace '+blind_string+freeze_string+freeze_r_string+' -t '+str(ntoys)+' -s '+str(seed) +' -n '+run_name

            if options.condor == True:
                tar_files = ['run_Stats.py',
                             'TwoDAlphabetClass.py',
                             'header.py',
                             'RpfHandler.py',
                             projDir+'/'+workspace_name]

                gof_toy_cmd = gof_toy_cmd.replace('GoodnessOfFit '+workspace_name,'GoodnessOfFit '+projDir+'/'+workspace_name)
                header.StatsForCondor(run_name,toy_dict,tar_files,[gof_toy_cmd])

            else:
                # header.executeCmd(gen_command,options.dryrun)
                header.executeCmd(gen_command)
                header.executeCmd(gof_toy_cmd)

        if not options.condor:
            if options.post:
                with header.cd('condor_'+run_name):
                    tmpdir = os.environ['CMSSW_BASE']+'/src/2DAlphabet/'+projDir+'/condor_'+run_name+'/tmp/'
                    header.executeCmd('mkdir '+tmpdir) 
                    if toyjobs > 1:
                        header.executeCmd('cat '+run_name+'*.tgz | tar zxvf - -i -C tmp/')
                        toyLimitTree = TChain('limit')
                        toyLimitTree.Add(tmpdir+'higgsCombine*.GoodnessOfFit.mH120.*.root') 
                    else:
                        header.executeCmd('tar xzvf *.tgz -C tmp/')
                        toyOutput = TFile.Open(tmpdir+'higgsCombine'+run_name+'.GoodnessOfFit.mH120.'+str(seed)+'.root')
                        toyLimitTree = toyOutput.Get('limit')
            else:
                toyOutput = TFile.Open('higgsCombine'+run_name+'.GoodnessOfFit.mH120.'+str(seed)+'.root')
                toyLimitTree = toyOutput.Get('limit')

            gStyle.SetOptStat(0)

            # Now to analyze the output
            # Get observation
            gofOutput = TFile.Open('higgsCombinegof_data.GoodnessOfFit.mH120.root')
            gofLimitTree = gofOutput.Get('limit')
            gofLimitTree.GetEntry(0)
            gofLimit = gofLimitTree.limit

            # Get toys
            # toyOutput = TFile.Open('higgsCombine'+run_name+'.GoodnessOfFit.mH120.'+str(seed)+'.root')
            # toyLimitTree = toyOutput.Get('limit')
            # toyLimitTree.Draw('>>limitEntryList','limit>1.0','entrylist')
            # limitEntryList = gDirectory.Get('limitEntryList')

            # toyLimitTree.SetEntryList(limitEntryList)
            # toyLimitTree.Print()
            # print toyLimitTree.GetMinimum('limit')
            # toyLimits = TH1F('hlimit','hlimit',25,800,1.2*max(toyLimitTree.GetMaximum('limit'),gofLimit))
            toyLimitTree.Draw('limit>>hlimit','limit>1.0 && limit<%s && limit != %s'%(gofLimit*2.0,gofLimit)) 
            toyLimits = gDirectory.Get('hlimit')
            time.sleep(1) # if you don't sleep the code moves too fast and won't perform the fit
            toyLimits.Fit("gaus")

            # Fit toys and derive p-value
            gaus = toyLimits.GetFunction("gaus")
            pvalue = 1-(1/gaus.Integral(-float("inf"),float("inf")))*gaus.Integral(-float("inf"),gofLimit)

            # Write out for reference
            out = open('gof_results.txt','w')
            out.write('Test statistic in data = '+str(gofLimit))
            out.write('Mean from toys = '+str(gaus.GetParameter(1)))
            out.write('Width from toys = '+str(gaus.GetParameter(2)))
            out.write('p-value = '+str(pvalue))

            # Extend the axis if needed
            if toyLimits.GetXaxis().GetXmax() < gofLimit:
                print 'Axis limit greater than GOF t value'
                binwidth = toyLimits.GetXaxis().GetBinWidth(1)
                xmin = toyLimits.GetXaxis().GetXmin()
                new_xmax = int(gofLimit*1.1)
                new_nbins = int((new_xmax-xmin)/binwidth)
                toyLimitTree.Draw('limit>>hlimitrebin('+str(new_nbins)+', '+str(xmin)+', '+str(new_xmax)+')','limit>0.001 && limit<1500') 
                toyLimits = gDirectory.Get('hlimitrebin')
                toyLimits.Fit("gaus")
                gaus = toyLimits.GetFunction("gaus")

            # Arrow for observed
            arrow = TArrow(gofLimit,0.25*toyLimits.GetMaximum(),gofLimit,0)
            arrow.SetLineWidth(2)

            # Legend
            leg = TLegend(0.1,0.7,0.4,0.9)
            leg.SetLineColor(kWhite)
            leg.SetLineWidth(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.AddEntry(toyLimits,"toy data","lep")
            leg.AddEntry(arrow,"observed = %.1f"%gofLimit,"l")
            leg.AddEntry(0,"p-value = %.2E"%decimal.Decimal(pvalue),"")

            # Draw
            cout = TCanvas('cout','cout',800,700)
            toyLimits.Draw('pez')
            arrow.Draw()
            leg.Draw()

            cout.Print('gof_plot.pdf','pdf')
            cout.Print('gof_plot.png','png')
            cout.SaveAs('gof_plot.root','root')

            if options.post:
                header.executeCmd('rm -r '+tmpdir)

    #########################################
    # Toy diagnositics and signal injection #
    #########################################
    if options.signalInjection != '':        
        if not options.plotOnly and not options.post:
            ###########################################
            # Now fit toys (send to condor if needed) #
            ###########################################
            # fit_command = 'combine -M FitDiagnostics'+mask_string.replace('=1','=0')+' -d initialFitWorkspace.root --snapshotName initialFit --toysFrequentist --skipBOnlyFit -t '+str(ntoys)+ ' --toysFile higgsCombine'+gen_name+'.GenerateOnly.mH120.'+str(seed)+'.root --rMin '+options.rMin+' --rMax '+options.rMax+' -n '+run_name
            fit_command = 'combine -M FitDiagnostics'+mask_string.replace('=1','=0')+' -d initialFitWorkspace.root --snapshotName initialFit --skipBOnlyFit --cminDefaultMinimizerStrategy 0 -t '+str(ntoys)+ ' --toysFile higgsCombine'+gen_name+'.GenerateOnly.mH120.'+str(seed)+'.root --rMin '+options.rMin+' --rMax '+options.rMax+' -n '+run_name+freeze_string
            if options.condor == True:
                tar_files = ['run_Stats.py',
                             'TwoDAlphabetClass.py',
                             'header.py',
                             'RpfHandler.py',
                             # projDir+'/'+workspace_name,
                             projDir+'/'+'initialFitWorkspace.root']

                this_gen_command = gen_command.replace('-d initialFitWorkspace.root','-d '+projDir+'/initialFitWorkspace.root')
                this_fit_command = fit_command.replace('-d initialFitWorkspace.root','-d '+projDir+'/initialFitWorkspace.root')

                header.StatsForCondor(run_name,toy_dict,tar_files,[this_gen_command,this_fit_command])

            else:
                header.executeCmd(gen_command,options.dryrun)
                header.executeCmd(fit_command,options.dryrun)

        ################################
        # Plot and save out the result #
        ################################
        if not options.condor:
            # If need to post-process a condor task
            if options.post:
                with header.cd('condor_'+run_name):
                    tmpdir = os.environ['CMSSW_BASE']+'/src/2DAlphabet/'+projDir+'/condor_'+run_name+'/tmp/'
                    header.executeCmd('mkdir '+tmpdir)
                    if toyjobs > 1:
                        header.executeCmd('cat '+run_name+'*.tgz | tar zxvf - -i -C tmp/')
                        tree_fit_sb = TChain('tree_fit_sb')
                        tree_fit_sb.Add(tmpdir+'fitDiagnostics'+run_name+'_*.root')
                    else:
                        header.executeCmd('tar xzvf *.tgz -C tmp/')
                        post_file = TFile.Open(tmpdir+'fitDiagnostics'+run_name+'.root')
                        tree_fit_sb = post_file.Get('tree_fit_sb')
            else:
                post_file = TFile.Open('fitDiagnostics'+run_name+'.root')
                tree_fit_sb = post_file.Get('tree_fit_sb')

            # Final plotting
            result_can = TCanvas('sigpull_can','sigpull_can',800,700)

            # if float(expectSignal) >= 0:
            tree_fit_sb.Draw("(r-"+expectSignal+")/(rHiErr*(r<"+expectSignal+")+rLoErr*(r>"+expectSignal+"))>>sigpull(20,-5,5)","fit_status>=0")
            tree_fit_sb.Draw("(r-"+expectSignal+")>>sigstrength(20,-1,1)","fit_status>=0")
            # else:
            #     tree_fit_sb.Draw("(r-1*("+expectSignal+"))/(rHiErr*(r<"+expectSignal+")+rLoErr*(r>"+expectSignal+"))>>sigpull(20,-5,5)","fit_status>=0")
            #     tree_fit_sb.Draw("(r-1*("+expectSignal+"))>>sigstrength(20,-1,1)","fit_status>=0")
            hsigpull = gDirectory.Get('sigpull')
            hsignstrength = gDirectory.Get('sigstrength')

            hsigpull.Fit("gaus","L")
            hsigpull.SetTitle(run_name)
            hsigpull.GetXaxis().SetTitle('(r-'+expectSignal+')/rErr')
            result_can.cd()
            hsigpull.Draw('pe')
            result_can.Print(run_name+'_sigpull.png','png')

            hsignstrength.Fit("gaus","L")
            hsignstrength.SetTitle(run_name)
            hsignstrength.GetXaxis().SetTitle('r-'+expectSignal)
            result_can.cd()
            hsignstrength.Draw('pe')
            result_can.Print(run_name+'_sigstrength.png','png')

            if options.post:
                header.executeCmd('rm -r '+tmpdir)
 
# Can run against other models or one model against itself to do a b-only vs s+b comparison
if options.biasStudy !='' or options.ftest:

    if options.ftest =='fitAlt' or options.ftest == 'pvalue':
        with header.cd(altDir):
            # Make a prefit workspace from the data card
            # if options.workspace == '':
            altworkspace_name = 'stats_workspace.root'
            t2w_cmd = 'text2workspace.py --channel-masks -b card_'+altcard_tag+'.txt -o '+altworkspace_name
            # if not (os.path.exists(workspace_name)):# and options.ftest == 'pvalue'):
            header.executeCmd(t2w_cmd,options.dryrun)

    ##############
    # Bias study #
    ##############
    if options.biasStudy:
        print 'Not working currently'
        # if not options.plotOnly and not options.post:
        #     ######################################################################
        #     # Go into the alternative model directory and run the toy fits there #
        #     ######################################################################
        #     fit_command = 'combine -M FitDiagnostics -d '+workspace_name+' -t '+str(ntoys)+ ' --toysFile '+projDir+'higgsCombine'+gen_name+'.GenerateOnly.mH120.'+str(seed)+'.root --rMin '+options.rMin+' --rMax '+options.rMax+' -n '+run_name
        #     if options.condor == True:
        #         this_gen_command = gen_command.replace('-d '+workspace_name,'-d '+projDir+'/'+workspace_name)
        #         this_fit_command = fit_command.replace('-d '+workspace_name,'-d '+altDir+'/'+workspace_name)
        #         tar_files = ['run_Stats.py',
        #                      'TwoDAlphabetClass.py',
        #                      'header.py',
        #                      'RpfHandler.py',
        #                      projDir+'/'+workspace_name,
        #                      altDir+'/'+altworkspace_name]

        #         # Need to cd into base dir to build the condor workspace
        #         with header.cd(projDir): header.StatsForCondor(run_name,toyjobs,tar_files,[this_gen_command,this_fit_command])

        #     else:
        #         with header.cd(projDir):
        #             header.executeCmd(gen_command,options.dryrun)
        #         with header.cd(altDir):
        #             this_fit_command = fit_command.replace(projDir+'higgsCombine',altDir_depth+'/'+projDir+'higgsCombine')
        #             header.executeCmd(this_fit_command,options.dryrun)

        # ############
        # # Plotting #
        # ############
        # if not options.condor:
        #     with header.cd(projDir):
        #         # If need to post-process a condor task
        #         if options.post:
        #             with header.cd(projDir+'/condor_'+run_name):
        #                 if toyjobs > 1:
        #                     header.executeCmd('cat '+run_name+'*.tgz | tar zxvf - -i')
        #                     tree_fit_sb = TChain('tree_fit_sb')
        #                     tree_fit_sb.Add('fitDiagnostics'+run_name+'_*.root')
        #                 else:
        #                     header.executeCmd('tar xzvf *.tgz')
        #                     post_file = TFile.Open('fitDiagnostics'+run_name+'.root')
        #                     tree_fit_sb = post_file.Get('tree_fit_sb')
        #         else:
        #             post_file = TFile.Open('fitDiagnostics'+run_name+'.root')
        #             tree_fit_sb = post_file.Get('tree_fit_sb')

        #     ################################
        #     # Plot and save out the result #
        #     ################################
        #     tree_fit_sb.Draw("(r-"+expectSignal+")/(rHiErr*(r<"+expectSignal+")+rLoErr*(r>"+expectSignal+"))>>sigpull(20,-5,5)","fit_status>=0")
        #     tree_fit_sb.Draw("(r-"+expectSignal+")>>sigstrength(20,-2,2)","fit_status>=0")

        #     hsigpull = gDirectory.Get('sigpull')
        #     hsignstrength = gDirectory.Get('sigstrength')

        #     hsigpull.Fit("gaus","L")
        #     hsigpull.SetTitle(run_name)
        #     hsigpull.GetXaxis().SetTitle('(r-'+expectSignal+')/rErr')
        #     result_can.cd()
        #     hsigpull.Draw('pe')
        #     result_can.Print(run_name+'_sigpull.png','png')

        #     hsignstrength.Fit("gaus","L")
        #     hsignstrength.SetTitle(run_name)
        #     hsignstrength.GetXaxis().SetTitle('r-'+expectSignal)
        #     result_can.cd()
        #     hsignstrength.Draw('pe')
        #     result_can.Print(run_name+'_sigstrength.png','png')

        #     header.executeCmd('rm condor_'+run_name+'/fitDiagnostics*.root condor_'+run_name+'/higgsCombine*.root')

    ##########
    # F test #
    ##########
    elif options.ftest:
        ftestseed = '123456'
        
        # Steps:
        # 1) Run GOF in base and alt models
        # 2) Run frequentist toy generation (from data) from base model
        # 3) Run GOF for both models over toys
        # 4) Compare results from steps 1 and 3
        base_fit_cmd = 'combine -d '+workspace_name+' -M GoodnessOfFit --algo saturated --rMax 10.0 --rMin -10.0 -n FTest '+blind_string+freeze_string+freeze_r_string
        if options.ftest == 'pvalue' or options.ftest == 'fitAlt': alt_fit_cmd = 'combine -d '+altworkspace_name+' -M GoodnessOfFit --algo saturated --rMax 10.0 --rMin -10.0 -n FTest '+blind_string+freeze_string+freeze_r_string

        if options.ftest == 'generate':#not options.plotOnly and options.ftest != 'post':
            # Base
            with header.cd(projDir):
                header.executeCmd(gen_command.replace(str(seed),ftestseed)) # Generation
            
        elif options.ftest == 'pvalue':
            with header.cd(projDir):
                header.executeCmd(base_fit_cmd) 
            with header.cd(altDir):
                header.executeCmd(alt_fit_cmd)

        elif options.ftest == 'fitMain':
            with header.cd(projDir):
                header.executeCmd(base_fit_cmd) 
                # Base fit to toys
                base_toy_fit_cmd = 'combine -d '+workspace_name+' -M GoodnessOfFit --rMax 10.0 --rMin -10.0 --algo saturated '+blind_string+freeze_string+freeze_r_string+' -t '+str(ntoys)+' --toysFile higgsCombineFTest.GenerateOnly.mH120.'+ftestseed+'.root -n FTestToyFits -s '+str(ftestseed)

                if options.condor:
                    base_toy_fit_cmd = base_toy_fit_cmd.replace('-d '+workspace_name,'-d '+projDir+'/'+workspace_name).replace('--toysFile higgs','--toysFile '+projDir+'/higgs')
                    tar_files = ['run_Stats.py',
                             'TwoDAlphabetClass.py',
                             'header.py',
                             'RpfHandler.py',
                             projDir+'/'+workspace_name,
                             projDir+'/higgsCombineFTest.GenerateOnly.mH120.'+ftestseed+'.root']
                             # altDir+'/'+altworkspace_name]

                    header.StatsForCondor(run_name,toy_dict,tar_files,[base_toy_fit_cmd,'ls'],files_to_grab=['higgsCombineFTestToyFits.GoodnessOfFit.mH120.123456.root'])

                else:
                    header.executeCmd(base_toy_fit_cmd)

        elif options.ftest == 'fitAlt':
            with header.cd(altDir):
                header.executeCmd(alt_fit_cmd)
                # Alt fit to toys
                alt_toy_fit_cmd = 'combine -d '+altworkspace_name+' -M GoodnessOfFit --rMax 10.0 --rMin -10.0 --algo saturated '+blind_string+freeze_string+freeze_r_string+' -t '+str(ntoys)+' --toysFile '+altDir_depth+projDir+'/higgsCombineFTest.GenerateOnly.mH120.'+ftestseed+'.root -n '+run_name+' -s '+str(ftestseed)

                if options.condor:
                    alt_toy_fit_cmd = alt_toy_fit_cmd.replace('-d '+workspace_name,'-d '+altDir+'/'+workspace_name).replace('--toysFile '+altDir_depth+projDir+'/higgs','--toysFile '+projDir+'/higgs')
                    tar_files = ['run_Stats.py',
                             'TwoDAlphabetClass.py',
                             'header.py',
                             'RpfHandler.py',
                             # projDir+'/'+workspace_name,
                             altDir+'/'+altworkspace_name,
                             projDir+'/higgsCombineFTest.GenerateOnly.mH120.'+ftestseed+'.root']

                    header.StatsForCondor(run_name,toy_dict,tar_files,[alt_toy_fit_cmd,'ls'],files_to_grab=['higgsCombine'+run_name+'.GoodnessOfFit.mH120.123456.root'])

                else:
                    header.executeCmd(alt_toy_fit_cmd)

        if options.ftest == 'post' or options.ftest == 'pvalue':
            # Do some basic checks before wasting compute time
            base_nrpf_params, base_nbins = header.ftestInfoLookup(header.projInfoLookup(projDir,card_tag))
            alt_nrpf_params, alt_nbins = header.ftestInfoLookup(header.projInfoLookup(altDir,altcard_tag))
            base_fit_filename = 'higgsCombineFTest.GoodnessOfFit.mH120.root'
            toy_fit_filename_alt = 'higgsCombine'+run_name+'.GoodnessOfFit.mH120.'+ftestseed+'.root'
            toy_fit_filename_main = 'higgsCombineFTestToyFits.GoodnessOfFit.mH120.'+ftestseed+'.root'

            # If ran with condor, go grab the output jobs
            if options.condor: 
                with header.cd(projDir+'/condor_'+run_name):
                    header.executeCmd('mkdir tmp/')
                    header.executeCmd('tar -xzvf '+run_name+'.tgz -C tmp/')
                header.executeCmd('mv '+projDir+'/condor_'+run_name+'/'+toy_fit_filename_main+' '+projDir+'/')

                with header.cd(altDir+'/condor_'+run_name):
                    header.executeCmd('mkdir tmp/')
                    header.executeCmd('tar -xzvf '+run_name+'.tgz -C tmp/')
                header.executeCmd('mv '+altDir+'/condor_'+run_name+'/'+toy_fit_filename_alt+' '+altDir+'/')

            # If the number of bins in the two models doesn't match, specify which to use or quit
            if base_nbins != alt_nbins:
                error_input = raw_input('ERROR: number of bins in the two models does not match (%i vs %i). Please repeat back which number you would like to calculate for or enter any other string to abort.'%(base_nbins,alt_nbins))
                if int(error_input) == base_nbins:
                    ftest_nbins = base_nbins
                elif int(error_input) == alt_nbins:
                    ftest_nbins = alt_nbins
                else:
                    print 'Quitting...'
                    quit()
            else:
                ftest_nbins = base_nbins


            base_fstat = header.FStatCalc(projDir+"/"+base_fit_filename, altDir+"/"+base_fit_filename, base_nrpf_params, alt_nrpf_params, ftest_nbins)
            if len(base_fstat) == 0: base_fstat = [0.0]

            ftest_p1 = min(base_nrpf_params,alt_nrpf_params)
            ftest_p2 = max(base_nrpf_params,alt_nrpf_params)
            ftest_nbins = base_nbins
            fdist = TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(10,1.3*base_fstat[0]))
            fdist.SetParameter(0,1)
            fdist.SetParameter(1,ftest_p2-ftest_p1)
            fdist.SetParameter(2,ftest_nbins-ftest_p2)

            if options.ftest == 'pvalue':
                pval = fdist.Integral(0.0,base_fstat[0])
                print 'P-value: %s'%pval
                pval_file = open('ftest_pval_'+projDir+'_vs_'+altDir+'.txt','w')
                pval_file.write(str(pval))
                pval_file.close()

                # Now we plot
                c = TCanvas('c','c',800,600)    
                c.SetLeftMargin(0.12) 
                c.SetBottomMargin(0.12)
                c.SetRightMargin(0.1)
                c.SetTopMargin(0.1)
                ftestHist_nbins = 30
                ftestHist = TH1F(tag+"Fhist",tag+"Fhist",ftestHist_nbins,0,max(10,1.3*base_fstat[0]))
                # ftestHist_cut = TH1F(tag+"Fhist_cut",tag+"Fhist cut",ftestHist_nbins,0,2*base_fstat[0])
                ftestHist.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
                ftestHist.GetXaxis().SetTitleSize(0.025)
                ftestHist.GetXaxis().SetTitleOffset(2)
                # # ftestHist.GetYaxis().SetTitle("Pseudodatasets")
                ftestHist.GetYaxis().SetTitleOffset(0.85)
                
                # for toyval in toys_fstat:
                #     ftestHist.Fill(toyval)
                #     if toyval > base_fstat[0]:
                #         ftestHist_cut.Fill(toyval)

                # ftestHist.SetMarkerStyle(20)
                ftestHist.Draw("pez")
                ftestobs  = TArrow(base_fstat[0],0.25,base_fstat[0],0)
                ftestobs.SetLineColor(kBlue+1)
                ftestobs.SetLineWidth(2)
                # ftestHist_cut.SetLineColor(kViolet-10)
                # ftestHist_cut.SetFillColor(kViolet-10)
                # ftestHist_cut.Draw("histsame")

                fdist.Draw('same')
                # int_fdist = TF1("int_fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(max(toys_fstat),base_fstat[0])+1)
                # int_fdist.SetParameter(0,1)
                # int_fdist.SetParameter(1,ftest_p2-ftest_p1)
                # int_fdist.SetParameter(2,ftest_nbins-ftest_p2)
                # pval = int_fdist.Integral(0.0,base_fstat[0])
              
                # ftestHist.Draw("pezsame")
                ftestobs.Draw()
                tLeg = TLegend(0.6,0.73,0.89,0.89)
                tLeg.SetLineColor(kWhite)
                tLeg.SetLineWidth(0)
                tLeg.SetFillStyle(0)
                tLeg.SetTextFont(42)
                tLeg.SetTextSize(0.03)
                # tLeg.AddEntry(ftestHist,"toy data","lep")
                tLeg.AddEntry(ftestobs,"observed = %.1f"%base_fstat[0],"l")
                #tLeg.AddEntry('',"p-value = %.2f"%(pval),"f")
                tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
                tLeg.Draw("same")

                # c.cd()
                model_info = TPaveText(0.2,0.6,0.4,0.8,"brNDC")
                model_info.AddText('p1 = '+projDir.split('pol')[1][0:4])
                model_info.AddText('p2 = '+altDir.split('pol')[1][0:4])
                model_info.AddText("p-value = %.2f"%(pval))
                model_info.Draw('same')
                
                latex = TLatex()
                latex.SetTextAlign(11)
                latex.SetTextSize(0.06)
                latex.SetTextFont(62)
                latex.SetNDC()
                latex.DrawLatex(0.12,0.91,"CMS")
                latex.SetTextSize(0.05)
                latex.SetTextFont(52)
                # if options.isData:
                latex.DrawLatex(0.65,0.91,"Preliminary")
                # else:
                #     l.DrawLatex(0.23,0.91,"Simulation")
                latex.SetTextFont(42)
                # latex.DrawLatex(0.76,0.91,"%.1f fb^{-1}"%options.lumi)
                latex.SetTextFont(52)
                latex.SetTextSize(0.045)
                # c.SaveAs(projDir+'/ftesst_vs_'+alttag+".root")
                c.SaveAs(projDir+'/ftest_vs_'+alttag+"_notoys.pdf")
                c.SaveAs(projDir+'/ftest_vs_'+alttag+"_notoys.png")

            elif options.ftest == 'post':
                print 'Analyzing F-test results...'
                toys_fstat = header.FStatCalc(projDir+"/"+toy_fit_filename_main, altDir+"/"+toy_fit_filename_alt, base_nrpf_params, alt_nrpf_params, ftest_nbins)

                toy_pass = 0
                for toyval in toys_fstat:
                    # print 'toys_fstat vs base_fstat:',toyval,base_fstat[0]
                    if base_fstat[0] > toyval:
                        toy_pass+=1

                # pval = 1
                # if len(toys_fstat) > 0:
                #     pval = float(toy_pass)/float(len(toys_fstat))
                #     print "F Test p-value",pval


                print "passing toys/number of toys = " + str(float(toy_pass)/float(len(toys_fstat)))

                # Now we plot
                c = TCanvas('c','c',800,600)    
                c.SetLeftMargin(0.12) 
                c.SetBottomMargin(0.12)
                c.SetRightMargin(0.1)
                c.SetTopMargin(0.1)
                ftestHist_nbins = 30
                ftestHist = TH1F(tag+"Fhist",tag+"Fhist",ftestHist_nbins,0,max(max(toys_fstat),base_fstat[0])+1)
                ftestHist_cut = TH1F(tag+"Fhist_cut",tag+"Fhist cut",ftestHist_nbins,0,max(max(toys_fstat),base_fstat[0])+1)
                ftestHist.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
                ftestHist.GetXaxis().SetTitleSize(0.025)
                ftestHist.GetXaxis().SetTitleOffset(2)
                ftestHist.GetYaxis().SetTitle("Pseudodatasets")
                ftestHist.GetYaxis().SetTitleOffset(0.85)
                
                for toyval in toys_fstat:
                    ftestHist.Fill(toyval)
                    if toyval > base_fstat[0]:
                        ftestHist_cut.Fill(toyval)

                ftestHist.SetMarkerStyle(20)
                ftestHist.Draw("pez")
                ftestobs  = TArrow(base_fstat[0],0.25*ftestHist.GetMaximum(),base_fstat[0],0)
                ftestobs.SetLineColor(kBlue+1)
                ftestobs.SetLineWidth(2)
                ftestHist_cut.SetLineColor(kViolet-10)
                ftestHist_cut.SetFillColor(kViolet-10)
                ftestHist_cut.Draw("histsame")

                fdist = TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(max(toys_fstat),base_fstat[0])+1)
                fdist.SetParameter(0,ftestHist.Integral()*((max(max(toys_fstat),base_fstat[0])+1)/ftestHist_nbins))
                fdist.SetParameter(1,ftest_p2-ftest_p1)
                fdist.SetParameter(2,ftest_nbins-ftest_p2)
                fdist.Draw('same')
                int_fdist = TF1("int_fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(max(toys_fstat),base_fstat[0])+1)
                int_fdist.SetParameter(0,1)
                int_fdist.SetParameter(1,ftest_p2-ftest_p1)
                int_fdist.SetParameter(2,ftest_nbins-ftest_p2)
                pval = int_fdist.Integral(0.0,base_fstat[0])
              
                ftestHist.Draw("pezsame")
                ftestobs.Draw()
                tLeg = TLegend(0.6,0.6,0.89,0.89)
                tLeg.SetLineColor(kWhite)
                tLeg.SetLineWidth(0)
                tLeg.SetFillStyle(0)
                tLeg.SetTextFont(42)
                tLeg.AddEntry(ftestHist,"toy data","lep")
                tLeg.AddEntry(ftestobs,"observed = %.1f"%base_fstat[0],"l")
                tLeg.AddEntry(ftestHist_cut,"p-value = %.2f"%(pval),"f")
                tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
                tLeg.Draw("same")

                # c.cd()
                model_info = TPaveText(0.2,0.6,0.4,0.8,"brNDC")
                model_info.AddText('p1 = '+projDir.split('pol')[1][0:4])
                model_info.AddText('p2 = '+altDir.split('pol')[1][0:4])
                model_info.Draw('same')
                
                latex = TLatex()
                latex.SetTextAlign(11)
                latex.SetTextSize(0.06)
                latex.SetTextFont(62)
                latex.SetNDC()
                latex.DrawLatex(0.12,0.91,"CMS")
                latex.SetTextSize(0.05)
                latex.SetTextFont(52)
                # if options.isData:
                latex.DrawLatex(0.23,0.91,"Preliminary")
                # else:
                #     l.DrawLatex(0.23,0.91,"Simulation")
                latex.SetTextFont(42)
                # latex.DrawLatex(0.76,0.91,"%.1f fb^{-1}"%options.lumi)
                latex.SetTextFont(52)
                latex.SetTextSize(0.045)
                # c.SaveAs(projDir+'/ftesst_vs_'+alttag+".root")
                c.SaveAs(projDir+'/ftest_vs_'+alttag+".pdf")
                c.SaveAs(projDir+'/ftest_vs_'+alttag+".png")
            if options.condor:
                header.executeCmd('rm -r '+altDir+'/condor_'+run_name+'/tmp/')
                header.executeCmd('rm -r '+projDir+'/condor_'+run_name+'/tmp/')   
