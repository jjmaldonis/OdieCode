
# Counts the atoms of each type in an xyz file
# You need specify in here on lines 15-17 which
# atomic numbers are in the xyz file

with open('asq_npt_heat_10000000.xyz') as f:
    lines = f.readlines()

one = 0
two = 0
three = 0
for i in range(2,len(lines)):
    lines[i] = lines[i].split()
    if(len(lines[i]) == 1): break
    if(lines[i][0].strip() == '40'): one = one + 1
    if(lines[i][0].strip() == '13'): two = two + 1
    if(lines[i][0].strip() == '29'): three = three + 1

print(one)
print(two)
print(three)

