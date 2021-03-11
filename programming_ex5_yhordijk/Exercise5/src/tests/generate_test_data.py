import random
import math
import os

folder = os.getcwd()



#generate test data for compress subroutine
def compress(a, m):
	indices = []
	for i in range(len(a)):
		if a[i] >= m:
			indices.append(i+1)
	return indices

modes = {
1: {'n': 30, 'm': 4, 'a': [1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,10,9,8,7,6,5,4,3,2,1]},
2: {'n': 10000, 'm': 22.4, 'a': [44.8*random.random() for i in range(10000)]},
3: {'n': 10000, 'm': 2, 'a': [2*random.random() for i in range(10000)]},
4: {'n': 10000, 'm': 4, 'a': [4 + random.random() for i in range(10000)]},
5: {'n': 10000, 'm': 4, 'a': [18 * random.random() for i in range(10000)]},
6: {'n': 10000, 'm': -1.5, 'a': [-3 * random.random() for i in range(10000)]},
7: {'n': 10000, 'm': 1, 'a': [(-1)**i for i in range(10000)]},
8: {'n': 10000, 'm': 0.6, 'a': [i*math.sin(2*math.pi*i/10000) for i in range(10000)]},
9: {'n': 10000, 'm': 0.6, 'a': [math.sin(100*math.pi*i/10000) for i in range(10000)]},
10: {'n': 10000, 'm': 50, 'a': [(i%(10000//100)) for i in range(10000)]},
11: {'n': 10000, 'm': 0.75, 'a': [8*((i/10000 - .5))**3 + (random.random()*2-1)/6 for i in range(10000)]}
}

for i, mode in modes.items():
	with open(folder + fr'\compress\test{("00", "0", "")[int(math.log10(i))]}{i}', 'w+') as f:
		n = mode['n']
		m = mode['m']
		a = mode['a']
		indices = compress(a,m)
		f.write(f'{n} {len(indices)} {m}\n')
		[f.write(f'{i}\n') for i in a]
		[f.write(f'{i}\n') for i in indices]





#generate test data for smooth subroutine
def smooth(a, m):
	smoothed = a.copy()
	for i in range(len(a)):
		if m < i + 1 < len(a)-m:
			s = 0
			for j in range(-m,m+1):
				s += a[i+j]
			smoothed[i] = s/(2*m+1)
	return smoothed

modes = {
1: {'n': 6, 'm': 1, 'a': [2,5,1,7,4,3]},
2: {'n': 10000, 'm': 5000, 'a': [random.random() for i in range(10000)]},
3: {'n': 1000, 'm': 1000, 'a': [random.random() for i in range(1000)]},
4: {'n': 500, 'm': 1000, 'a': [random.random() for i in range(500)]},
5: {'n': 10000, 'm': 5000, 'a': [(-1)**i for i in range(10000)]},
6: {'n': 10000, 'm': 5000, 'a': [i*math.sin(2*math.pi*i/10000) for i in range(10000)]},
7: {'n': 10000, 'm': 5000, 'a': [math.sin(100*math.pi*i/10000) for i in range(10000)]},
8: {'n': 10000, 'm': 5000, 'a': [(i%(10000//100)) for i in range(10000)]},
9: {'n': 10000, 'm': 5000, 'a': [(0,5000)[i%5000==(5000//2)] for i in range(10000)]},
10: {'n': 10000, 'm': 1000, 'a': [8*((i/10000 - .5))**3 + (random.random()*2-1)/6 for i in range(10000)]}
}

for i, mode in modes.items():
	with open(folder + fr'\smooth\test{("00", "0", "")[int(math.log10(i))]}{i}', 'w+') as f:
		n = mode['n']
		m = mode['m']
		a = mode['a']
		f.write(f'{n} {m}\n')
		[f.write(f'{i} {k}\n') for i, k in zip(a, smooth(a, m))]



#generate test data for filter subroutine
def filter(a, m):
	filtered = a.copy()
	for i in range(len(a)):
		filtered[i] = min(m, max(-m, a[i]))
	return filtered

modes = {
1: {'n': 6, 'm': 3, 'a': [2,-5,1,7,-4,3]},
2: {'n': 10000, 'm': 0.5, 'a': [random.random() for i in range(10000)]},
3: {'n': 10000, 'm': 0, 'a': [random.random()-0.5 for i in range(10000)]},
4: {'n': 10000, 'm': 0.25, 'a': [random.random()-0.5 for i in range(10000)]},
5: {'n': 10000, 'm': 0.5, 'a': [random.random()+0.5 for i in range(10000)]},
}

for i, mode in modes.items():
	with open(folder + fr'\filter\test{("00", "0", "")[int(math.log10(i))]}{i}', 'w+') as f:
		n = mode['n']
		m = mode['m']
		a = mode['a']
		f.write(f'{n} {m}\n')
		[f.write(f'{i} {k}\n') for i, k in zip(a, filter(a, m))]

