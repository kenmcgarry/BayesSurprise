s = "ABACBABAABBCBABA"
d={}
MINCNT=2
MINLEN=2

for sublen in range(MINLEN,int(len(s)/MINCNT)):
    for i in range(0,len(s)-sublen):
        sub = s[i:i+sublen]
        cnt = s.count(sub)
        if cnt >= MINCNT and sub not in d:
             d[sub] = cnt
             
             
