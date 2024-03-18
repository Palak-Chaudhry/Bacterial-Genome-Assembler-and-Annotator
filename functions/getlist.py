import os
entries = os.listdir('./raw_fq/')
entries = [x.split('.')[0] for x in entries]
entries = [x.split('_')[0] for x in entries]
entries = set(entries)

with open('list.txt', 'w') as file:
    for entry in entries:
        file.write(entry + '\n')

with open('listx.txt', 'w') as file:
    for entry in entries:
        file.write(entry + ' ')