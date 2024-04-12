using Plots

begin
	Nₓ = 20
	Lₓ = 1.0
	δx = Lₓ/Nₓ 
	xs = δx/2:δx:Lₓ

	Nₜ = 1000
	t = 1.0
	δt = t/Nₜ

	U = 1

	D = 0.05

	T₀ = sin.(2π * xs)
end

function advection_diffusion(T_0, δt, δx, U)
	N = length(T_0)
	T = similar(T_0)

	for i in 1:N
        i_prev = i == 1 ? N : i - 1
        i_next = i == N ? 1 : i + 1

        adv = -(U * δt) * (T_0[i_next] - T_0[i_prev]) / (2 * δx)
        diff = D * (δt / (δx)^2) * (T_0[i_next] + T_0[i_prev] - 2 * T_0[i])
        T[i] = T_0[i] + adv + diff
	end

	return T
end

f(x,t) = sin(2 * π * (x - U * t)) * exp( - 4* π^2 * D *t)

function Tₑ(x,t)
	T = []
	for xᵢ ∈ x
		push!(T, f(xᵢ,t))
	end
	return T
end

typeof(xs)

T = [T₀]

for i ∈ 1:Nₜ
	push!(T, advection_diffusion(last(T), δt, δx, U))
end

@gif for i ∈ 1:3:Nₜ
	plot(xs, T[i], ylims = [-1,1])
	plot!(xs, Tₑ(xs, i * δt), ylims = [-1,1])
end
