import ROOT
from ROOT import *

import subprocess
from optparse import OptionParser

def generate(name,distEq,xVar,yVar,nevents):

    # Gaussian
    if distEq.find('gauss') != -1:
        meanx = RooConstVar(name+'_meanx',name+'_meanx',15)
        sigmax = RooConstVar(name+'_sigmax',name+'_sigmax',2)

        meany = RooConstVar(name+'_meany',name+'_meany',10)
        sigmay = RooConstVar(name+'_sigmay',name+'_sigmay',2)


        if distEq == 'gauss':
            dummyPDFx = RooGaussian(name+'x',name+'x',xVar,meanx,sigmax)
            dummyPDFy = RooGaussian(name+'y',name+'y',yVar,meany,sigmay)

        elif distEq == 'gaussUp':
            tailx = RooConstVar('tailx','tailx',0.5)
            taily = RooConstVar('taily','taily',0.5)
            dummyPDFx = RooNovosibirsk(name+'x',name+'x',xVar,meanx,sigmax,tailx)
            dummyPDFy = RooNovosibirsk(name+'y',name+'y',yVar,meany,sigmay,taily)

        elif distEq == 'gaussDown':
            tailx = RooConstVar('tailx','tailx',-0.5)
            taily = RooConstVar('taily','taily',-0.5)
            dummyPDFx = RooNovosibirsk(name+'x',name+'x',xVar,meanx,sigmax,tailx)
            dummyPDFy = RooNovosibirsk(name+'y',name+'y',yVar,meany,sigmay,taily)

        dummyPDF = RooProdPdf(name+'_PDF',name+'_PDF',RooArgList(dummyPDFx,dummyPDFy))

    # Generic
    else:
        dummyPDF = RooGenericPdf(name+'_PDF',distEq,RooArgList(xVar,yVar))

    dummyRDS = dummyPDF.generate(RooArgSet(xVar,yVar),nevents)
    dummyRDH = RooDataHist(name+'_RDH',name+'_RDH',RooArgSet(xVar,yVar),dummyRDS)
    dummyTH2 = dummyRDS.createHistogram(xVar,yVar,24,20,'',name)

    dummyTH2.SetName(name)

    return dummyTH2


if __name__ == '__main__':
    gStyle.SetOptStat(0)

    #########################################################
    #                       Options                         #
    #########################################################
    parser = OptionParser()
    # Input and what to run
    parser.add_option('-N', '--name', type='string', action='store',
                    default   =   '',
                    dest      =   'name',
                    help      =   'Name for output')
    parser.add_option('-F', '--failform', type='string', action='store',
                    default   =   '',
                    dest      =   'failform',
                    help      =   'Equation for failing distribution (in x and y)')
    parser.add_option('-P', '--passform', type='string', action='store',
                    default   =   '',
                    dest      =   'passform',
                    help      =   'Equation for passing distribution (in x and y)')
    parser.add_option('-f', '--nevents_fail', type='int', action='store',
                    default   =   10000000,
                    dest      =   'nevents_fail',
                    help      =   'Number of events to generate in fail')
    parser.add_option('-p', '--nevents_pass', type='int', action='store',
                    default   =   10000000,
                    dest      =   'nevents_pass',
                    help      =   'Number of events to generate in pass')
    
    # parser.add_option('-s', '--neventsSig', type='integer', action='store',
    #                 default   =   1000,
    #                 dest      =   'neventsSig',
    #                 help      =   'Number of events to generate for signal')

    (options, args) = parser.parse_args()

    outfile = TFile('distributions/'+options.name+'.root','RECREATE')

    # Establish our axis variables
    xVar = RooRealVar('x','x',0,24)
    yVar = RooRealVar('y','y',0,20)

    if len(options.failform.split(',')) != len(options.passform.split(',')):
        print 'Number of fail forms does not match number of pass forms. Quiting'
        quit()
    for i,fail_form in enumerate(options.failform.split(',')):
        pass_form = options.passform.split(',')[i]
        if i == 0:
            thisname = 'data_obs'
        else:
            thisname = 'bkg'+str(i)
        bkg_fail = generate(thisname+'_fail',fail_form,xVar,yVar,options.nevents_fail)
        bkg_pass = generate(thisname+'_pass',pass_form,xVar,yVar,options.nevents_pass)
        bkg_fail.Write()
        bkg_pass.Write()

    signal_fail = generate('signal_fail','gauss',xVar,yVar,options.nevents_fail/100)
    signalUp_fail = generate('signalUp_fail','gaussUp',xVar,yVar,options.nevents_fail/100)
    signalDown_fail = generate('signalDown_fail','gaussDown',xVar,yVar,options.nevents_fail/100)
    signal_pass = generate('signal_pass','gauss',xVar,yVar,options.nevents_pass/20)
    signalUp_pass = generate('signalUp_pass','gaussUp',xVar,yVar,options.nevents_pass/20)
    signalDown_pass = generate('signalDown_pass','gaussDown',xVar,yVar,options.nevents_pass/20)

    signal_fail.Write()
    signal_pass.Write()

    signalUp_fail.Write()
    signalUp_pass.Write()

    signalDown_fail.Write()
    signalDown_pass.Write()

    outfile.Close()


