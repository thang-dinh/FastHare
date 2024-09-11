import fasthare as m

assert m.__version__ == '1.0.4'

ski1 = [(0, 1, 3)]
ski2 = [(0, 1, -3)]
# exp1.net
rh, map, sign, time = m.fasthare_reduction(ski1)
assert (len(rh) == 0) and (map[0] + map[1] == 0) and (sign[0] + sign[1] == 2) 

#exp2.net
rh, map, sign, time = m.fasthare_reduction(ski2)
assert (len(rh) == 0) and (map[0] + map[1] == 0) and (sign[0] + sign[1] == 0) 

#exp3.net
rh3  = []
map3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 
sign3= [-1, 1, 1, -1, -1, -1, -1, -1, -1, -1]
rh, map, sign, time = m.fasthare_reduction(file ="exp3.net",alpha=0.2, log_file="exp3.log")
assert (rh == rh3) and (map == map3) and (sign == sign3) 
print("All tests passed!")
