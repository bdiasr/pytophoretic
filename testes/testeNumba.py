import time
import numpy as np
from numba import njit

array = np.random.randint(1000, 10000, 10000000)

def preco_com_imposto(lista):
    aumento_total = 0
    for item in lista:
        if item > 5000:
            valor = item * 1.1
        else:
            valor = item
        aumento_total += valor - item
    return aumento_total

inicio = time.time()
print("Aumento:", preco_com_imposto(array))
fim = time.time()
print(fim - inicio)


@njit
def preco_com_imposto(lista):
    aumento_total = 0
    for item in lista:
        if item > 5000:
            valor = item * 1.1
        else:
            valor = item
        aumento_total += valor - item
    return aumento_total

inicio = time.time()
print("Aumento:", preco_com_imposto(array))
fim = time.time()
print(fim - inicio)