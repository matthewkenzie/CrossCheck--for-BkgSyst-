import os
import sys
import fnmatch

fitName = ['2pol','3pol','4pol','1exp','2exp','3exp','1pow','2pow','3pow','2lau','4lau','6lau','title']

for fit in fitName:
  htF = open('plots/Graphs/'+fit+'.html','w')
  htF.write('Author: Matthew Kenzie <br> \n')
  htF.write('<script language="Javascript"> \n document.write("Last Modified: " + document.lastModified + ""); \n </script> <br> \n ')
  htF.write('<center>\n <p> \n <font size="5">Background Model Studies</font>\n </p>\n')
  htF.write('<a href=anRes.log>Log file</a><br>\n')
  htF.write('<font size="3">Results for '+fit+' fit</font><br> \n')
  htF.write('Fit: \n')
  for gen in fitName:
    if gen=='title':
      continue
    htF.write('<a href=\"'+gen+'.html\">'+gen+'</a>\n')
  htF.write('<br>\n')
  for root,dirs,files in os.walk('plots/Graphs'):
    for filename in fnmatch.filter(files,'*.png'):
      if 'f'+fit in filename:
        htF.write('<a href='+filename+'><img height=\"400\" src=\"'+filename+'\"></a>\n')
  
