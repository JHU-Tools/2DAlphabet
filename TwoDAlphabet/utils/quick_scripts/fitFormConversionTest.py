RFVform = '(@1+@2*x+@3*x**2)'
shift = 3
TF1form = ''
lookingForDigit = False
for index,char in enumerate(RFVform):
    # print str(index) + ' ' + char
    if char == '@':
        TF1form+='['
        lookingForDigit = True
    elif lookingForDigit:
        if char.isdigit():
            if RFVform[index+1].isdigit() == False:     # if this is the last digit
                TF1form+=str(int(char)+shift)
        else:
            TF1form+=']'+char
            lookingForDigit = False
    else:
        TF1form+=char

print TF1form