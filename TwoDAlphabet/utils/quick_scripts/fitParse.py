fitForm = '5*(@0+(@1*(x))+(@2*x**2))(@3+@4*y)'   # MUST BE FACTORIZED LIKE (X PART)(YPART)

# Look for all opening and closing parentheses
openIndexes = []
closedIndexes = []
for index,char in enumerate(fitForm):
    if char == '(': # start looking for ")" after this "("
        openIndexes.append(index)
    if char == ')':
        closedIndexes.append(index)


# Now pair them by looking at the first in closedIndexes and finding the closest opening one (to the left)
# Remove the pair from the index lists and repeat
innerContent = []
for iclose in closedIndexes:
    diff = len(fitForm)     # max length to start because we want to minimize
    for iopen in openIndexes:
        if iclose > iopen:
            this_diff = iclose - iopen
            if this_diff < diff:
                diff = this_diff
                candidateOpen = iopen
    openIndexes.remove(candidateOpen)
    innerContent.append(fitForm[iclose-diff+1:iclose])


outerContent = []
for c in innerContent:
    keep_c = True
    for d in innerContent:
        if '('+c+')' in d and c != d:
            keep_c = False
            break
    if keep_c:
        outerContent.append(c)

if len(outerContent) != 2:
    print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
    quit()
else:
    for c in outerContent:
        if 'x' in c and 'y' not in c:
            xPart = c
        elif 'x' not in c and 'y' in c:
            yPart = c
        else:
            print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
            quit()

print xPart
print yPart
