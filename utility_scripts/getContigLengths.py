#!usr/bin/env python

import sys 

file_open = open(sys.argv[1],'r').readlines() 
file_out = open(sys.argv[2],'w') 

switch = 0 
for lines in file_open:
    if '>' in lines:
        if switch == 0:
            file_out.write(lines.strip() + "\t")
            seq_length = 0
            switch = 1
        elif switch == 1:
            file_out.write(str(seq_length) + "\n" + lines.strip() + "\t")
            seq_length = 0
    else:
        seq_length+=len(lines.strip())
file_out.write(str(seq_length))
