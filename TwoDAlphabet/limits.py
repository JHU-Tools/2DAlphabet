def MakeSmoothGraph(h2,h3):
    h2 = ROOT.TGraph(h2)
    h3 = ROOT.TGraph(h3)
    npoints = h3.GetN()
    h3.Set(2*npoints+2)
    for b in range(npoints+2):
        x1, y1 = (ROOT.Double(), ROOT.Double())
        if b == 0:
            h3.GetPoint(npoints-1, x1, y1)
        elif b == 1:
            h2.GetPoint(npoints-b, x1, y1)
        else:
            h2.GetPoint(npoints-b+1, x1, y1)
        h3.SetPoint(npoints+b, x1, y1)
    return h3

def Inter(g1,g2):
    xaxisrange = g1.GetXaxis().GetXmax()-g1.GetXaxis().GetXmin()
    xaxismin = g1.GetXaxis().GetXmin()
    xinters = []
    yinters = []
    for x in range(0,10000):
        xpoint = xaxismin + (float(x)/1000.0)*xaxisrange
        xpoint1 = xaxismin + (float(x+1)/1000.0)*xaxisrange
        Pr1 = g1.Eval(xpoint)
        Pr2 = g2.Eval(xpoint)
        Po1 = g1.Eval(xpoint1)
        Po2 = g2.Eval(xpoint1)
        if (Pr1-Pr2)*(Po1-Po2)<0:
            xinters.append(0.5*(xpoint+xpoint1))
            yinters.append(0.5*(Po1+Po2))
        
    if len(xinters) == 0:
        xinters = [-1]
        yinters = [-1]
    return xinters[0],yinters[0]