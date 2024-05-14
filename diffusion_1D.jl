using Plots

begin
  N = 21
  L = 1.0
  h = L / (N - 1)
  x = 0:h:L

  Nₜ = 1000
  t = 1.0
  τ = t / Nₜ
  T = τ:τ:t

  D = 0.05

  CFL = τ * D / h^2

  u₀ = zeros(N)
  u₀[1] = 1/h
end

G(x, t) = (1 / (sqrt(4 * D * π * t))) * exp.(-(x) .^ 2 ./ (4 * D * t))

function is_stable(CFL)
  return CFL < 1
end

function diffusion_explicit(u₀, CFL)
  N = length(u₀)
  u = similar(u₀)

  u[1] = u₀[1] + 2 * CFL * (u₀[2] - u₀[1]) # idea: use u₀[i-1] = u₀[i+1]

  for i in 2:N-1
    u[i] = u₀[i] + CFL * (u₀[i+1] + u₀[i-1] - 2 * u₀[i])
  end

  u[N] = u₀[N] + CFL * (u₀[N-1] - u₀[N]) # idea: use u₀[N+1] = u₀[N]

  return u
end

function diffusion_implicit(u₀, CFL)
  N = length(u₀)
  u = similar(u₀)

  A = zeros(N, N)
  b = zeros(N)

  A[1, 1] = 1 + CFL
  A[1, 2] = - CFL

  b[1] = (1 - CFL) * u₀[1] + CFL * u₀[2] 

  for i in 2:N-1
    A[i, i] = 1 + CFL
    A[i, i-1] = -0.5 * CFL
    A[i, i+1] = -0.5 * CFL
    b[i] = 0.5 * CFL * u₀[i+1] + (1 - CFL) * u₀[i] + 0.5 * CFL * u₀[i-1]
  end

  A[N, N] = 1 + 0.5 * CFL
  A[N, N-1] = -0.5 * CFL

  b[N] = (1 - 0.5 * CFL) * u₀[N] + 0.5 * CFL * u₀[N-1]

  u = A \ b

  return u
end

if is_stable(CFL)
  Uₑ = [u₀]
  Uᵢ = [u₀]

  for i ∈ 1:Nₜ
    push!(Uₑ, diffusion_explicit(last(Uₑ), CFL))
    push!(Uᵢ, diffusion_implicit(last(Uᵢ), CFL))
  end

  @gif for i ∈ 1:2:Nₜ
    plot(x, Uₑ[i], ylims=[0, 1 / h])
    plot!(x, Uᵢ[i])
    plot!(x, G(x, T[i]))
  end
else
  println("UNSTABLE")
end