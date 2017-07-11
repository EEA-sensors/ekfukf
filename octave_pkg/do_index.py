

with open ('../Contents.m', 'r') as f:
  index = f.read().strip()

import re
m = re.search(r'(?P<header>.*)%\s$(?P<content>.*)', index, re.S|re.I|re.M)

header = m.group('header')
index  = m.group('content').strip()

end_re   = re.compile (r'^%\s*$', re.M)
index = end_re.split(index)

topic_re = re.compile (r'^%\s([A-Z][a-z\s]+[a-zA-Z\s]+)$', re.M)
func_re  = re.compile (r'^%\s+([A-Z_0-9]+)\s+',re.M)

with open ('INDEX','w') as outf:
  outf.write('ekfukf >> Karmal filters\n')
  for section in index:
    if section:
      topic = topic_re.findall (section)
      if topic:
        outf.write('{}\n'.format(topic[0]))
        funcs = func_re.findall (section)
        [outf.write('  {}\n'.format(fn.lower())) for fn in funcs]
