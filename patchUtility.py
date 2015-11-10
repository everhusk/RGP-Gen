import sys

original = sys.argv[1]
dest     = sys.argv[2]

newcode = """elseif (hvar.eq.'ptpn') then
        v=ptpn(icomp)
        if (l) ptpn(icomp)=dbl
      elseif (hvar.eq.'txeos') then"""

with open(original, 'r') as utility_file:
  utility = utility_file.read();
  
utility = utility.replace("elseif (hvar.eq.'txeos') then",newcode)

with open(dest, 'w') as utility_file:
  utility_file.write(utility)

print('patch complete')
