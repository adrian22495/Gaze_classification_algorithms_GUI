import csv
#Eyetribe
f = open("Test1_1135582_201604121106014314.txt", "r")
f.readline()
f.readline()
f.readline()

csvsalida = open('Test1_1135582_201604121106014314.csv', 'w', newline='')
salida = csv.writer(csvsalida, delimiter=',')
for line in f:
    row = line.split(', ')
	
    #7 timestamp, 13 & 14 avg both eyes, 17 & 18 left eye, 23 & 24 right eye
    salida.writerow([row[7], row[13], row[14]])
        
del f, salida, row
csvsalida.close()
del csvsalida