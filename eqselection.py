

from mhps.earthquake import get_earthquake_list, get_total_excitations
import numpy as np
import os


excitation_file = 'Excitations.csv'

total_eq = get_total_excitations(excitation_file)
earthquake_generator = get_earthquake_list(excitation_file)

idx = int(input('Start index for storing : '))
print(total_eq)

for i in range(total_eq):
    ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)
    xg = xg.reshape(len(xg), 1)
    yg = yg.reshape(len(yg), 1)
    print('Time Step = %8.4f'%(dt*ndiv))
    print('EQ Duration = %8.4f'%(dt*ndiv*ndt))
    print('PGA in X = %8.4f g'%(max(abs(xg/9.81))))
    print('PGA in Y = %8.4f g'%(max(abs(yg/9.81))))
    t = np.arange(0, dt*ndt, dt )
    t = t.reshape(len(t),1)
    print(t.shape, xg.shape)
    file = os.path.join('EQs', 'EQ-' + str(idx) + '.csv')
    eq = np.hstack((t, xg/9.81, yg/9.81))
    np.savetxt(file, eq, delimiter=",")
    idx = idx + 1
    a = input('\nPress Enter to continue...\n')