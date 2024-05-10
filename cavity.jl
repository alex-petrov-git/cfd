using Plots
#using PlotlyJS

function Solve()
  Nx = 17
  Ny = 17
  MaxStep = 200
  ε = 0.1
  dt = 0.005
  t = 0.0
  h = 1.0 / (Nx - 1)
  MaxIt = 100
  Beta = 1.5
  MaxErr = 0.001

  # Parameters for SOR
  ψ = zeros(Nx, Ny)
  ψ₀ = zeros(Nx, Ny)
  ω = zeros(Nx, Ny)
  ω₀ = zeros(Nx, Ny)
  u = zeros(Nx, Ny)
  v = zeros(Nx, Ny)
  x = 0.0:h:1.0
  y = 0.0:h:1.0

  for istep in 1:MaxStep
    for iter in 1:MaxIt
      ψ₀ .= ψ
      for i in 2:Nx-1, j in 2:Ny-1
        ψ[i, j] = 0.25 * Beta * (ψ[i+1, j] + ψ[i-1, j] + ψ[i, j+1] + ψ[i, j-1] + h * h * ω[i, j]) + (1.0 - Beta) * ψ[i, j]
      end
      Err = 0.0
      for i in 1:Nx, j in 1:Ny
        Err += abs(ψ₀[i, j] - ψ[i, j])
      end
      if Err <= MaxErr
        break
      end
    end

    for i in 2:Nx-1, j in 2:Ny-1
      u[i, j] = (ψ[i, j+1] - ψ[i, j-1])/(2*h)
      v[i, j] = - (ψ[i+1, j] - ψ[i-1, j]) / (2 * h)
    end

    ω[2:Nx-1, 1] .= -2.0 * ψ[2:Nx-1, 2] / (h * h)  # vorticity on bdrys
    ω[2:Nx-1, Ny] .= -2.0 * ψ[2:Nx-1, Ny-1] / (h * h) .- 2.0 / h  # top wall
    ω[1, 2:Ny-1] .= -2.0 * ψ[2, 2:Ny-1] / (h * h)  # right wall
    ω[Nx, 2:Ny-1] .= -2.0 * ψ[Nx-1, 2:Ny-1] / (h * h)  # left wall
    ω₀ .= ω
    for i in 2:Nx-1, j in 2:Ny-1
      ω[i, j] += dt * (-0.25 * ((ψ[i, j+1] - ψ[i, j-1]) * (ω₀[i+1, j] - ω₀[i-1, j]) - (ψ[i+1, j] - ψ[i-1, j]) * (ω₀[i, j+1] - ω₀[i, j-1])) / (h * h) + ε * (ω₀[i+1, j] + ω₀[i-1, j] + ω₀[i, j+1] + ω₀[i, j-1] - 4.0 * ω₀[i, j]) / (h^2))
    end
    t += dt
  end

  return [x, y, ω, ψ, u, v]
end

x, y, ω, ψ, u, v = Solve()

# Создание контурного графика и заполнение

contour_plot = contourf(x, y, ω, alpha=0.5, color=:viridis, xlabel="X", ylabel="Y", title="Каверна (функция тока)")

#****************************************************************

#use PlotlyJS 

#surf1 = surface(x=x, y=y, z=ω, colorscale=:blues)
#surf2 = surface(x=x, y=y, z=ψ, colorscale=:blues)

#layout1 = Layout(scene=attr(xaxis_title="x", yaxis_title="y", zaxis_title="ω", title="Каверна (завихренность)"))
#layout2 = Layout(scene=attr(xaxis_title="x", yaxis_title="y", zaxis_title="ψ", title="Каверна (функция тока)"))

#plot(surf1, layout1)
#plot(surf2, layout2)

#****************************************************************

#c1 = Plots.heatmap(x, y, ω)
#c2 = Plots.heatmap(x, y, ψ)

#Plots.plot(c1, c2)