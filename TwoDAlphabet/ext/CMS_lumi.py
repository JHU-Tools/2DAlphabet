import ROOT as rt

# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#

cmsText     = "CMS";
cmsTextFont   = 61  

writeExtraText = True
extraText   = "Preliminary"
extraTextFont = 52 

lumiTextSize     = 0.35
lumiTextOffset   = 0.2

cmsTextSize      = 0.4
cmsTextOffset    = 0.15

relPosX    = 0.07
relPosY    = -0.05
relExtraDX = 2.8
relExtraDY = 0.3

extraOverCmsTextSize  = 0.76

lumi_13TeV = "137 fb^{-1}"
lumi_2016_13TeV = "35.9 fb^{-1}"
lumi_2017_13TeV = "41.5 fb^{-1}"
lumi_2018_13TeV = "60 fb^{-1}"
lumi_8TeV  = "19.7 fb^{-1}" 
lumi_7TeV  = "5.1 fb^{-1}"
lumi_sqrtS = ""

drawLogo      = False

def CMS_lumi(pad,  iPeriod=4,  iPosX=11, sim=False ):
    outOfFrame    = False
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    # if( iPeriod==1 ):
    #     lumiText += lumi_7TeV
    #     lumiText += " (7 TeV)"
    # elif ( iPeriod==2 ):
    #     lumiText += lumi_8TeV
    #     lumiText += " (8 TeV)"

    # elif( iPeriod==3 ):      
    #     lumiText = lumi_8TeV 
    #     lumiText += " (8 TeV)"
    #     lumiText += " + "
    #     lumiText += lumi_7TeV
    #     lumiText += " (7 TeV)"
    if ( iPeriod==1 ):
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==16 ):
        lumiText += lumi_2016_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==17 ):
        lumiText += lumi_2017_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==18 ):
        lumiText += lumi_2018_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==2 ):
        lumiText += lumi_2016_13TeV+' + '+lumi_2017_13TeV + ' + ' + lumi_2018_13TeV
        lumiText += " (13 TeV)"
    # elif ( iPeriod==7 ):
    #     if( outOfFrame ):lumiText += "#scale[0.85]{"
    #     lumiText += lumi_13TeV 
    #     lumiText += " (13 TeV)"
    #     lumiText += " + "
    #     lumiText += lumi_8TeV 
    #     lumiText += " (8 TeV)"
    #     lumiText += " + "
    #     lumiText += lumi_7TeV
    #     lumiText += " (7 TeV)"
    #     if( outOfFrame): lumiText += "}"
    # elif ( iPeriod==12 ):
    #     lumiText += "8 TeV"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
            
    # print lumiText

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    
    
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    

    latex.DrawLatex(1-r,1-t+cmsTextOffset*t,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+cmsTextOffset*t,cmsText)
  
    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo =  rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()          
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                if not sim: latex.DrawLatex(posX_+ relExtraDX*cmsTextSize*t, posY_-relExtraDY*cmsTextSize*t, extraText)
                else: latex.DrawLatex(posX_+relExtraDX*cmsTextSize*t, posY_-relExtraDY*cmsTextSize*t, extraText + ' simulation')
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        if not sim: latex.DrawLatex(posX_, posY_, extraText)
        else: latex.DrawLatex(posX_, posY_, extraText + ' simulation')

    pad.Update()

