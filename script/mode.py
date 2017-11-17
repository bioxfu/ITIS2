from collections import Counter

def mode(num):
	count_dict = dict(Counter(num).most_common())
	max_count = Counter(num).most_common(1)[0][1]
	modes = []
	for x in num:
		if count_dict[x] == max_count: 
			modes.append(x)
	mode = sorted(modes)[0]
	return(mode)

print mode([1087607, 1087628])
print mode([1220051,1220059,1220063,1220068])

