'''Script to use 2D Alphabet output to calculate the pull and saturated test statistic
per-bin over the full 2D space.'''

import sys,os,ROOT,math,array
import header

def stitchHists(name,thisHistList,blinded=[]):
    # Required that thisHistList be in order of desired stitching
    # `blinded` is a list of the index of regions you wish to skip/blind
    # Version of what's in header.py. Difference is that this deduces the xbinning
    # from the input histograms.
    xbins = []
    for h in thisHistList:
        for ix in range(1,h.GetNbinsX()+1):
            xbins.append(h.GetXaxis().GetBinLowEdge(ix))
    xbins.append(thisHistList[-1].GetXaxis().GetBinUpEdge(thisHistList[-1].GetNbinsX()))

    ybins = []
    for iy in range(1,thisHistList[0].GetNbinsY()+1):
        ybins.append(thisHistList[0].GetYaxis().GetBinLowEdge(iy))
    ybins.append(thisHistList[0].GetYaxis().GetBinUpEdge(thisHistList[0].GetNbinsY()))

    axbins = array.array('d',xbins)
    aybins = array.array('d',ybins)
    stitched_hist = ROOT.TH2F(name,name,len(xbins)-1,axbins,len(ybins)-1,aybins)

    bin_jump = 0
    for i,h in enumerate(thisHistList):
        if i in blinded:
            bin_jump += thisHistList[i].GetNbinsX()
            continue
        
        for ybin in range(1,h.GetNbinsY()+1):
            for xbin in range(1,h.GetNbinsX()+1):
                stitched_xindex = xbin + bin_jump

                stitched_hist.SetBinContent(stitched_xindex,ybin,h.GetBinContent(xbin,ybin))
                stitched_hist.SetBinError(stitched_xindex,ybin,h.GetBinError(xbin,ybin))

        bin_jump += thisHistList[i].GetNbinsX()

    return stitched_hist

def getSaturated(h_data,h_bkg):
    h_saturated = h_data.Clone()
    h_saturated.Reset()
    h_saturated.SetName(h_saturated.GetName().replace('_data','_saturated'))
    
    for ix in range(1,h_saturated.GetNbinsX()+1):
        for iy in range(1,h_saturated.GetNbinsY()+1):
            f = h_bkg.GetBinContent(ix,iy)
            d = h_data.GetBinContent(ix,iy)

            if f > 0 and d > 0 :
                s = f - d + d * math.log(d/f)
            else:
                s = 0

            h_saturated.SetBinContent(ix,iy,s)

    name_parts = h_saturated.GetName().split('_')
    new_title = '%s %s %s'%(name_parts[0],name_parts[1],'saturated')
    h_saturated.SetTitle(new_title)

    return h_saturated

def getPull(h_data,h_bkg):
    h_pull = h_data.Clone()
    h_pull.Reset()
    h_pull.SetName(h_pull.GetName().replace('_data','_pull'))
    
    for ix in range(1,h_pull.GetNbinsX()+1):
        for iy in range(1,h_pull.GetNbinsY()+1):
            f = h_bkg.GetBinContent(ix,iy)
            d = h_data.GetBinContent(ix,iy)
            ferr = h_bkg.GetBinError(ix,iy)
            derr = h_data.GetBinError(ix,iy)

            if ferr > 0 or derr > 0:
                p = (d-f)/math.sqrt(ferr**2+derr**2)
            else:
                p = 0
            
            h_pull.SetBinContent(ix,iy,p)

    name_parts = h_pull.GetName().split('_')
    new_title = '%s %s %s'%(name_parts[0],name_parts[1],'pull')
    h_pull.SetTitle(new_title)

    return h_pull


if __name__ == "__main__":
    projDir = sys.argv[1]
    subdirs = [x for x in os.listdir(projDir) if os.path.isdir(projDir+'/'+x)]

    postfitshapes = ROOT.TFile.Open(projDir+'/postfitshapes_b.root')

    to_plot_saturated = []
    titles_saturated = []
    to_plot_pull = []
    titles_pull = []

    sat_integral = 0

    for r in ['pass','fail']:
        for sd in subdirs:
            data_to_stitch = []
            bkg_to_stitch = []
            for cat in ['LOW','SIG','HIGH']:
                data_to_stitch.append(postfitshapes.Get('%s_%s_%s_postfit/data_obs'%(r,cat,sd)))
                bkg_to_stitch.append(postfitshapes.Get('%s_%s_%s_postfit/TotalBkg'%(r,cat,sd)))

            h_data = stitchHists('%s_%s_data'%(r,sd),data_to_stitch,blinded=[1])
            h_bkg = stitchHists('%s_%s_bkg'%(r,sd),bkg_to_stitch,blinded=[1])

            if 'saturated' in sys.argv:
                sat = getSaturated(h_data,h_bkg)
                to_plot_saturated.append(sat)
                titles_saturated.append(sat.GetTitle())
                sat_integral += sat.Integral()
            if 'pull' in sys.argv:
                pull = getPull(h_data,h_bkg)
                to_plot_pull.append(pull)
                titles_pull.append(pull.GetTitle())

    print ('Integral of saturated plots: %s'%(sat_integral))

    if len(to_plot_saturated) > 0: header.makeCan('twoD_saturated',projDir,to_plot_saturated,xtitle='m_{h} (GeV)',ytitle='m_{hh} (GeV)',titles=titles_saturated,datastyle='lego',ztitle='f-d+d*ln(d/f)')
    if len(to_plot_pull) > 0: header.makeCan('twoD_pull',projDir,to_plot_pull,xtitle='m_{h} (GeV)',ytitle='m_{hh} (GeV)',titles=titles_pull,datastyle='colz',ztitle='(data-bkg)/#sigma')