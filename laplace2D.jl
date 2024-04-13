using PlotlyJS

function laplace2D_Gauss_Seilder(nx::Int, ny::Int, max_iter::Int, tolerance::Float64)
  # Длины сторон прямоугольника
  Lx = 1.0
  Ly = 1.0

  # Шаги сетки
  dx = Lx / nx
  dy = Ly / ny

  a = 0.5 * dx^2/(dx^2 + dy^2)
  b = 0.5 * dy^2/(dx^2 + dy^2)

  # Создание начального приближения
  phi = zeros(nx, ny + 2)  # Создаем массив для потенциала (с учетом граничных узлов)

  # Установка граничных условий
  phi[:, 1] .= 0  # Нижняя граница

  # Граничное условие справа
  phi[:, end] .= 1

  # Численное решение уравнения Лапласа
  for iter in 1:max_iter
    # Создание копии текущего приближения для проверки сходимости
    phi_old = copy(phi)

    # Цикл по внутренним узлам сетки (исключая граничные узлы)
    for i in 2:nx-1, j in 2:ny+1
      # Аппроксимация вторых производных оператором Лапласа
      phi[i, j] = b * (phi[i+1, j] + phi[i-1, j]) + a * (phi[i, j+1] + phi[i, j-1])
    end

    # Проверка сходимости
    max_residual = maximum(abs.(phi - phi_old))
    if max_residual < tolerance
      println("Gauss_Seilder: Сходимость достигнута на шаге $iter")
      break
    end
  end

  return phi[2:end-1, 2:end-1]
end

function laplace2D_SOR(nx::Int, ny::Int, max_iter::Int, tolerance::Float64, β::Float64)
  # Длины сторон прямоугольника
  Lx = 1.0
  Ly = 1.0

  # Шаги сетки
  dx = Lx / nx
  dy = Ly / ny

  a = 0.5 * dx^2 / (dx^2 + dy^2)
  b = 0.5 * dy^2 / (dx^2 + dy^2)

  # Создание начального приближения
  phi = zeros(nx, ny + 2)  # Создаем массив для потенциала (с учетом граничных узлов)

  # Установка граничных условий
  phi[:, 1] .= 0  # Нижняя граница

  # Граничное условие справа
  phi[:, end] .= 1

  # Численное решение уравнения Лапласа
  for iter in 1:max_iter
    # Создание копии текущего приближения для проверки сходимости
    phi_old = copy(phi)

    # Цикл по внутренним узлам сетки (исключая граничные узлы)
    for i in 2:nx-1, j in 2:ny+1
      # Аппроксимация вторых производных оператором Лапласа
      phi[i, j] = (b * (phi[i+1, j] + phi[i-1, j]) + a * (phi[i, j+1] + phi[i, j-1])) * β + (1 - β) * phi[i, j]
    end

    # Проверка сходимости
    max_residual = maximum(abs.(phi - phi_old))
    if max_residual < tolerance
      println("SOR: Сходимость достигнута на шаге $iter")
      break
    end
  end

  return phi[2:end-1, 2:end-1]
end

function plot_solution(phi; nx, ny)
  x = range(0, stop=1, length=nx)
  y = range(0, stop=1, length=ny)

  surf = surface(x=x, y=y, z=phi, colorscale=:blues, name="Численное решение")
  
  layout = Layout(scene=attr(xaxis_title="x", yaxis_title="y", zaxis_title="Phi", title="Решение уравнения Лапласа"))

  plot(surf, layout)
end
# Параметры для численного решения
nx = 50  # Число узлов по оси x
ny = 50  # Число узлов по оси y
max_iter = 10000  # Максимальное число итераций
tolerance = 1e-6  # Порог сходимости
β = 1.9

# Вызов функции для решения уравнения Лапласа
phi_Gauss_Seilder = laplace2D_Gauss_Seilder(nx, ny, max_iter, tolerance)
phi_SOR = laplace2D_SOR(nx, ny, max_iter, tolerance, β)

plot_solution(phi_SOR; nx=nx, ny=ny)

