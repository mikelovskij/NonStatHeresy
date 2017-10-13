def tabmaker(tabcontent, theader = False, lheader = False):
    tstr = "<table id=\"t01\">\n"
    for (row, i) in zip(tabcontent, xrange(len(tabcontent))):
        tstr +='<tr>\n'
        if theader & (i == 0):
            for elem in row:
                tstr +='<th>\n'
                tstr += elem + '\n'
                tstr += '</th>\n'
        else:
            for (elem, j) in zip(row, xrange(len(row))):
                if lheader & (j == 0):
                    tstr += '<th>\n'
                    tstr += elem + '\n'
                    tstr += '</th>\n'
                else:
                    tstr += '<td>\n'
                    tstr += elem + '\n'
                    tstr += '</td>\n'
        tstr += '</tr>\n'
    tstr += '</table>\n'
    return tstr
