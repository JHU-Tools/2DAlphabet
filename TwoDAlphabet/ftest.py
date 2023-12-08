from ROOT import TFile

def FtestInfoLookup(projInfoDict):
    nrpfparams = 0
    nbins = 0
    for k in projInfoDict.keys():
        this_nrpfparams = len(projInfoDict[k]['rpfVarNames'])
        this_nbins = (len(projInfoDict[k]['full_x_bins'])-1) * (len(projInfoDict[k]['newYbins'])-1)
        if projInfoDict[k]['blindedFit'] == True:
            this_nbins = this_nbins - ( (len(projInfoDict[k]['newXbins']['SIG'])-1) * (len(projInfoDict[k]['newYbins'])-1) )

        nrpfparams+=this_nrpfparams
        nbins+=this_nbins
    
    return nrpfparams,nbins

def FstatCalc(filename1,filename2,p1,p2,n):
    print ('Calculating F statistic')
    # Flip flop to make sure p2 is always greater than p1 (more parameters should always fit better)
    if p1 > p2:
        p1, p2 = p2, p1
        filename1, filename2 = filename2, filename1
    print ('Files: ',filename1,filename2)
    print ('Parameters: p1 %f, p2 %f, n %f'%(p1,p2,n))

    # Get limit trees from each file
    file1 = TFile.Open(filename1)
    tree1 = file1.Get("limit")
    file2 = TFile.Open(filename2)
    tree2 = file2.Get("limit")
    diffs=[]
    # print 'Entries ',tree1.GetEntries(),tree2.GetEntries()

    # Loop over entries and calculate the F statistics
    for i in range(0,tree1.GetEntries()):
        tree1.GetEntry(i)
        tree2.GetEntry(i)
        # print 'Limits ',tree1.limit,tree2.limit
        if tree1.limit-tree2.limit>0:
            F = (tree1.limit-tree2.limit)/(p2-p1)/(tree2.limit/(n-p2))
            # print 'Entry ',i, ":", tree2.limit, "-", tree1.limit, "=", tree2.limit-tree1.limit, "F =", F
            # if F < 50: 
            diffs.append(F)
        else:
            print ('WARNING in calculation of F statistic for entry %i. limit1-limit2 <=0 (%f - %f)' %(i,tree1.limit,tree2.limit))
            diffs.append(0)
    # print 'Diffs F stat: ',diffs
    return diffs
