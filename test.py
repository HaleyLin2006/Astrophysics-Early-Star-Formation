import json

with open('convert.txt') as json_file:
    json_decoded = json.load(json_file)

print([float(i) for i in json_decoded["0,1,2,3,4,5,6,7,8,9,10"].split(",")])
json_decoded[','.join([str(i) for i in range(11)])] = ','.join([str(i) for i in range(5)])

with open('convert.txt', 'w') as json_file:
    json.dump(json_decoded, json_file)
