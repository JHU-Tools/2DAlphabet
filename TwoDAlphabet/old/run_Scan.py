'''Run a likelihood scan along a specified parameter (default r)
and plot the output.
'''
import ROOT, header, array, os

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--projDir", dest="projDir",default='',
                help="Home of the project - has the cards, fit results, etc")
parser.add_option("-t", dest="toys",default='',
                help="Number of toys to do")
parser.add_option("-p", dest="parameter",default='r',
                help="Parameter to scan. Defaults to r.")
parser.add_option("--toysFile", dest="toysFile",default='1',
                help="File holding toys")
parser.add_option("-n", dest="name",default='',
                help="Extra name identifier.")
parser.add_option("--seed", dest="seed",default='',
                help="Seed - fed to combine")
parser.add_option('--rMin', metavar='F', type='string', action='store',
                default =   '-5',
                dest    =   'rMin',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option('--rMax', metavar='F', type='string', action='store',
                default =   '5',
                dest    =   'rMax',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option("--skipScan", 
                action="store_true", dest="skipScan", default=False,
                help="Skip scan if just plotting")
(options, args) = parser.parse_args()

projDir = options.projDir # home of the workspace - has the cards, fit results, etc
tag = projDir.split('/')[0]
runname = 'Scan%s'%options.name
rrange = [int(options.rMin),int(options.rMax)]

if options.toys == '':
    raise ValueError('Toy option empty. Scanning data fit not currently supported.')

if os.path.exists(options.toysFile):
    toysFile = options.toysFile
elif os.path.exists(projDir+'/'+options.toysFile):
    toysFile = projDir+'/'+options.toysFile
else:
    raise Exception('Cannot find file %s'%options.toysFile)

with header.cd(projDir):
    if not options.skipScan: header.executeCmd('combine -M MultiDimFit -d card_%s.txt --algo grid -P %s --floatOtherPOIs=0 --setParameterRanges %s=%s,%s --cminDefaultMinimizerStrategy 0 --saveNLL -n %s -t %s --toysFile %s'%(tag,options.parameter,options.parameter,rrange[0],rrange[1],runname,options.toys,options.toysFile))

    header.executeCmd('mkdir %s_plots/'%runname)
    f_out = ROOT.TFile.Open('%s_plots/all.root'%runname,'RECREATE')
    f = ROOT.TFile.Open('higgsCombine%s.MultiDimFit.mH120.123456.root'%runname)
    for i in range(1,int(options.toys)+1):
        TLimit = f.Get('limit')
        
        r = []
        nll = []
        for ientry in range(1,TLimit.GetEntries()):
            TLimit.GetEntry(ientry)
            if TLimit.deltaNLL*2 !=0 and TLimit.iToy == i:
                nll.append(2*TLimit.deltaNLL)
                r.append(TLimit.r)
                print ('r = %s, nll = %s'%(TLimit.r,TLimit.deltaNLL))

        r = array.array('d',r)
        nll = array.array('d',nll)
        graph = ROOT.TGraph(len(r),r,nll)
        graph.SetName('toy_scan_%s'%i)
        graph.SetTitle('Toy %s'%i)

        c = ROOT.TCanvas('c','c',800,700)
        graph.Draw('AC')
        c.Print('%s_plots/toy_scan_%s.png'%(runname,i),'png')

        f_out.cd()
        graph.Write()

        del graph
        del c

    f_out.Close()
