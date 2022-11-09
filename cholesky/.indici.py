import numpy as np

n = 5 #shape matrice quadrata

external = 0 #variabile per contare quanti cicli esterni devo fare
internal = 2 * n - 1 - 1

aux = 0
for col in range(2 * n - 1):
    if col < n-1:
        row = 0
        to_print = str(row) + ' ' +  str(col)
        aux += 1
        for z in range(1, int(np.floor(col/2))+1):
            to_print += ' | ' + str(row+z) + ' ' +  str(col-z)
        
        print(to_print)

    else:
        row = n-1
        to_print = str(col-aux) + ' ' +  str(row)
        for z in range(1, int(np.floor(internal/2))+1):
            to_print += ' | ' + str(col-aux+z) + ' ' +  str(row-z)

        print(to_print)

    internal -= 1   
    external += 1

print(f'cicli esterni: {external}')


print()



n = 6 #shape matrice quadrata

external = 0 #variabile per contare quanti cicli esterni devo fare
internal = 2 * n - 1 - 1

aux = 0
for row in range(2 * n - 1):
    if row < n-1:
        col = 0
        to_print = str(row) + ' ' +  str(col)
        aux += 1
        for z in range(1, int(np.floor(row/2))+1):
            to_print += ' | ' + str(row-z) + ' ' +  str(col+z)
        
        print(to_print)

    else:
        col = n-1
        to_print = str(col) + ' ' +  str(row-aux)
        for z in range(1, int(np.floor(internal/2))+1):
            to_print += ' | ' + str(col-z) + ' ' +  str(row-aux+z)

        print(to_print)

    internal -= 1   
    external += 1

print(f'cicli esterni: {external}')