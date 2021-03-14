import random
import math
import os
import numpy as np
import numpy.lib.recfunctions as rec

folder = os.getcwd()


#define cases here with given rows, cols and parameter c and of course array a
cases = {
1: {'rows': 6, 'cols': 3, 'c':1, 'a': np.round(np.random.random(size=(6,3)),5)},
2: {'rows': 6, 'cols': 3, 'c':2, 'a': np.round(np.random.random(size=(6,3)),5)},
3: {'rows': 6, 'cols': 3, 'c':3, 'a': np.round(np.random.random(size=(6,3)),5)},
4: {'rows': 12, 'cols': 3, 'c':1, 'a': np.round(np.random.randint(0,4,size=(12,3)),5)},
5: {'rows': 12, 'cols': 3, 'c':1, 'a': np.round(np.random.random(size=(12,3)),5)},
6: {'rows': 12, 'cols': 3, 'c':1, 'a': np.round(np.random.random(size=(12,3)),5)},
7: {'rows': 12, 'cols': 3, 'c':1, 'a': np.round(np.random.random(size=(12,3)),5)},
}

#set column 1 of the 5th test-case to all zero
cases[5]['a'][:,0] = 0.0
#sort column 1 of 6th test-case
cases[6]['a'] = cases[6]['a'][np.argsort(cases[6]['a'][:,cases[6]['c']-1], kind='stable')]
#sort column 1 of 7th test-case and flip it
cases[7]['a'] = cases[7]['a'][np.argsort(cases[7]['a'][:,cases[7]['c']-1], kind='stable')][::-1,:]

#loop over cases
for i, case in cases.items():
	#create new file
	with open(folder + fr'\test{("00", "0", "")[int(math.log10(i))]}{i}', 'w+') as f:
		#write number of rows and cols and c for fortran to know array size
		f.write(f"{case['rows']} {case['cols']} {case['c']}\n")
		#write array to txt, append newline char after each row
		[[f.write((f'{el}\t') + '\n' * ((j+1) == case['cols'])) for j, el in enumerate(row)] for row in case['a']]
		#sort array and write to txt
		#first argsort column c - 1 and then use those indices to rearrange original array
		case['a'] = case['a'][np.argsort(case['a'][:,case['c']-1], kind='stable')]
		[[f.write((f'{el}\t') + '\n' * ((j+1) == case['cols'])) for j, el in enumerate(row)] for row in case['a']]