import numpy as np

n = 5 #shape matrice quadrata

external = 0 #variabile per contare quanti cicli esterni devo fare
internal = 2 * n - 1 - 1

aux = 0
for i in range(2 * n - 1):
    if i < n-1:
        j = 0
        to_print = str(j) + ' ' +  str(i)
        aux += 1
        for z in range(1, int(np.floor(i/2))+1):
            to_print += ' | ' + str(j+z) + ' ' +  str(i-z)
        
        print(to_print)

    else:
        j = n-1
        to_print = str(i-aux) + ' ' +  str(j)
        for z in range(1, int(np.floor(internal/2))+1):
            to_print += ' | ' + str(i-aux+z) + ' ' +  str(j-z)

        print(to_print)

    internal -= 1   
    external += 1

print(f'cicli esterni: {external}')