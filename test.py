import json

details = {"id": 5, 
           "name": "haley", 
           "class": 8, 
           "email": "haley@pynative.com"}

with open('convert.txt', 'a') as convert_file:
	convert_file.write(json.dumps(details))

# reading the data from the file
studentsList = []
print("Started Reading JSON file which contains multiple JSON document")
with open('convert.txt') as f:
    for jsonObj in f:
        studentDict = json.loads(jsonObj)
        studentsList.append(studentDict)

print("Printing each JSON Decoded Object")
for student in studentsList:
    print(type(student))
    print(student["id"], student["name"], student["class"], student["email"])