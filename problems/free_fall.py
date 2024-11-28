import math
# Данные задачи
m = 23
g = 9.8
F_resistance = 13 # сила сопротивления при падении на 700 метров
H = 700 
v0 = 0  
x0 = 0  

Err = 0.00001 # желаемая ошибка в силе
Err_current = 1

# solver parameters
dt = 0.0001
a = 1
V = [v0]

while Err_current > Err:
    h = 0
    V = [v0]
    while h < H:
        V.append(V[-1] + dt*(g - (a/m)*(V[-1])**2))
        h += dt*(V[-1] + V[-2])/2
    Err_current = abs(F_resistance - a*(V[-1])**2) 
    a = a * (F_resistance / (a * (V[-1])**2))

print(math.sqrt(m*g/a))