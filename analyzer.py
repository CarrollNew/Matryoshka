with open('hmmscan_clear') as file2:
	f2 = file2.read().split('\n')
trdict = dict()
for i2 in range(len(f2)):
	if f2[i2]:
		trdict[f2[i2].split(' ')[16]] = None
coords = 10000, 12000
with open('headers') as file:
	f = file.read().split('\n')
Borders1, Borders2 = False, False

for i in range(len(f)):
	if f[i]:
		id = f[i].split(' ')[0][1:]
		c1 = f[i].split(' ')[2]
		c2 = f[i].split(' ')[4]
		if id in trdict:
			print(id, c1, c2, "Transposase")
		else:
			print(id, c1, c2, "Non-transposase")
