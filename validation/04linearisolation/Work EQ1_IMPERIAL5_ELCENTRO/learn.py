


import csv


line = []
line.append('FUNDAMENTAL SUPERSTRUCTURE TIME PERIOD: ' + str(2.0) +'SEC')
line.append('RATIO OF BASE MASS TO FLOOR:' + str(1.0))
line.append('PERIOD OF ISOLATED STRUCTURE')
line.append('MASS OF BASE RAFT')
line.append('DAMPING RATIO OF SLIDING SYSTEM')



with open('people.csv', 'a') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow(line)
    writer.writerow('\n')