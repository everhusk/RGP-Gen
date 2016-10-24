from __future__ import print_function
import sys

original = sys.argv[1]
dest     = sys.argv[2]

newcode = """elseif (hvar.eq.'ptpn') then
        v=ptpn(icomp)
        if (l) ptpn(icomp)=dbl
      elseif (hvar.eq.'txeos') then"""

def addPTPN():
  with open(original, 'r') as utility_file:
    utility = utility_file.read()
    
  # don't patch if already patched
  # (unlikely as we get the original file)
  if "elseif (hvar.eq.'ptpn') then" not in utility:
    utility = utility.replace("elseif (hvar.eq.'txeos') then",newcode)

  with open(dest, 'w') as utility_file:
    utility_file.write(utility)

addPTPN()
print('REFPROP: Patch complete', end="")