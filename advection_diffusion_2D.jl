using Plots

function advection_diffusion(T_0, δt, δl, U, V, D)
    Nx, Ny = size(T_0)
    T = similar(T_0)

    for i in 1:Nx, j in 1:Ny
        # Вычисляем индексы соседних ячеек с учетом граничных условий
        i_prev = i == 1 ? Nx : i - 1
        i_next = i == Nx ? 1 : i + 1
        j_prev = j == 1 ? Ny : j - 1
        j_next = j == Ny ? 1 : j + 1

        # Вычисляем значения переноса и диффузии вдоль обеих осей
        adv_x = U * (T_0[i_next, j] - T_0[i_prev, j]) / (2 * δl)
        adv_y = V * (T_0[i, j_next] - T_0[i, j_prev]) / (2 * δl)
        diff_x = D * (T_0[i_next, j] + T_0[i_prev, j] - 2 * T_0[i, j]) / (δl^2)
        diff_y = D * (T_0[i, j_next] + T_0[i, j_prev] - 2 * T_0[i, j]) / (δl^2)

        # Обновляем значение температуры в текущей ячейке
        T[i, j] = T_0[i, j] - δt * (adv_x + adv_y) + δt * (diff_x + diff_y)
    end

    return T
end

function initial_temperature(xs, ys, x0, y0, σ)
    A = 1.0  # Амплитуда
    T₀ = zeros(length(xs), length(ys))

    for i in 1:length(xs)
        for j in 1:length(ys)
            T₀[i, j] = A * exp(-((xs[i] - x0)^2 + (ys[j] - y0)^2) / (2σ^2))
        end
    end

    return T₀
end

begin
	N = 100
	L = 1.0
	δl = L/N
	xs = δl/2:δl:L
	ys = δl/2:δl:L

	Nₜ = 3000
	t = 1.0
	δt = t/Nₜ

	U = 1

	V = -1

	D = 0.05

	T₀ = initial_temperature(xs, ys, 0.5, 0.5, 0.1)
end

function is_stable(U, V, D, δt, δl)
    for L in 0:π/50:2*π
        for M in 0:π/50:2*π
            Re_q = 1 - 4 * D * δt / δl^2 * (sin(L/2)^2 + sin(M/2)^2)
            Im_q = - (U * δt / (2 * δl)) * sin(L) - (V * δt / (2 * δl)) * sin(M)
            mod_q = sqrt(Re_q^2 + Im_q^2)
            if mod_q > 1
                return false
            end
        end
    end
    return true
end

T = [T₀]

if is_stable(U, V, D, δt, δl)
	for i ∈ 1:Nₜ
		push!(T, advection_diffusion(last(T), δt, δl, U, V, D))
	end
else 
	println("UNSTABLE")
end

function plot_animation(T, xs, ys, δl)
    anim = @animate for i in 1:10:length(T)
        heatmap(xs, ys, T[i], xlabel="x", ylabel="y", color=:viridis, c=:blues, aspect_ratio=:equal)
    end

    return anim
end

animation = plot_animation(T, xs, ys, δl)

gif(animation, "temperature_animation.gif", fps = 150)
