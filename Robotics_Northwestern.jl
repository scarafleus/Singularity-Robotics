### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 36844d82-3f80-11eb-244a-83efb39d115d
begin
	using Plots
	using PlutoUI
	using LinearAlgebra
end

# ╔═╡ 54276a22-3f80-11eb-2538-e18635928ecf
md"# Robotics

`scarafleus`, `blueSky`, `Kirite`, `Albert Wigmore`. December 2020 - ?"

# ╔═╡ 350a9f42-3f85-11eb-06ae-2556cc9aa4fd
md"###### Axis"

# ╔═╡ fb3fedce-4d55-11eb-2879-b755a16cde4e
md"Commonly used vectors in their homogeneous space representation."

# ╔═╡ 327e9936-3f85-11eb-1bf9-d1a2def49b41
begin
	xhat::Array{Float64, 1} = [1.0, 0.0, 0.0, 1.0]
	yhat::Array{Float64, 1} = [0.0, 1.0, 0.0, 1.0]
	zhat::Array{Float64, 1} = [0.0, 0.0, 1.0, 1.0]
	zero::Array{Float64, 1} = [0.0, 0.0, 0.0, 1.0]
end

# ╔═╡ c7f5bac4-4da1-11eb-333f-63859c65f8f7
md"#### Matrices for computations"

# ╔═╡ 83fd510c-42ab-11eb-347c-4f9625869adc
md"###### Special Matrices"

# ╔═╡ dcdfa036-4d55-11eb-1026-1763c6107e31
md"Computes the identity matrix of dimension `n` by `n`."

# ╔═╡ 8d4bdea4-42ab-11eb-1f3a-e5cea01bcf3f
eye(n) = Matrix(1.0I, n, n)

# ╔═╡ f2a1fbbc-4d55-11eb-0ff1-41e28c57ad71
md"Computes the identity matrix of dimension `n` by `m`."

# ╔═╡ faa7d0d0-4d36-11eb-38e8-371dc6654639
eye(n, m) = Matrix(1.0I, n, m)

# ╔═╡ f0952168-3f82-11eb-37f2-b17628f1044c
md"###### Rotation Matrix"

# ╔═╡ 9be3e8c6-4d55-11eb-1949-43d88b18b0db
md"Computes a rotation matrix by combining a rotation axis `ω` with an angle `θ`."

# ╔═╡ cd374f48-3f82-11eb-26f6-8dfbc7216940
function rotationMatrix(ω, θ)
	ω = ω/sqrt(ω[1]^2 + ω[2]^2 + ω[3]^2)
	return rotationMatrix(ω[1], ω[2], ω[3], θ)
end

# ╔═╡ 6e303e66-4d55-11eb-3c6f-59cfcbf9211c
md"Computes a rotation matrix by combining a rotation axis given by `ω1`, `ω2`, `ω3` with an angle `θ`."

# ╔═╡ 5e89ba96-3f81-11eb-0a35-e9d3928427b8
function rotationMatrix(ω1, ω2, ω3, θ)
	cθ = cos(θ)
	sθ = sin(θ)
	
	R = zeros(3, 3)
	R[1, 1] = cθ + ω1^2*(1 - cθ)
	R[2, 1] = ω1*ω2*(1 - cθ) + ω3*sθ
	R[3, 1] = ω1*ω3*(1 - cθ) - ω2*sθ
	R[1, 2] = ω1*ω2*(1 - cθ) - ω3*sθ
	R[2, 2] = cθ + ω2^2*(1 - cθ)
	R[3, 2] = ω2*ω3*(1 - cθ) + ω1*sθ
	R[1, 3] = ω1*ω3*(1 - cθ) + ω2*sθ
	R[2, 3] = ω2*ω3*(1 - cθ) - ω1*sθ
	R[3, 3] = cθ + ω3^2*(1 - cθ)
	return R
end

# ╔═╡ 30385498-3f83-11eb-1584-792e07e24dc7
rotationMatrix(zhat, π/6)

# ╔═╡ 0829b69a-3f83-11eb-3f65-0f5ab0b35aa9
md"###### Transformation Matrix"

# ╔═╡ 176e3a1a-4d55-11eb-0dd2-27d5cdac3b20
md"Computes the 4x4 dimensional transformation matrix from combining a rotation matrix `R` with an offset vector `p`."

# ╔═╡ 152049ea-3f83-11eb-18dd-8d667de8625d
function transformationMatrix(R, p)
	T = zeros(4, 4)
	T[1:3, 1:3] = R
	T[1:3, 4] = p[1:3]
	T[4, 4] = 1
	return T
end

# ╔═╡ 3abe0ebe-4d55-11eb-05f8-35ce4ffefd66
md"Computes the 4x4 dimensional transformation matrix from combining a rotation matrix given by a rotation axis `ω` and an angle `θ` with an offset vector `p`."

# ╔═╡ 09009752-3f86-11eb-1e79-3f70a93d7575
transformationMatrix(ω, θ, p) = transformationMatrix(rotationMatrix(ω, θ), p)

# ╔═╡ c2826700-4d55-11eb-19cf-7b814acfca64
md"Extracts the rotation matrix from a transformation matrix `T`."

# ╔═╡ 0d673ae0-42ac-11eb-05ec-7798c0839217
R(T) = T[1:3, 1:3]

# ╔═╡ cf770c18-4d55-11eb-362f-a138a8b6c353
md"Extracts the offset vector from a transformation matrix `T`."

# ╔═╡ 1db602d0-42ac-11eb-0992-735dc410de36
p(T) = T[1:3, 4]

# ╔═╡ dee42ada-3f83-11eb-2c39-83a2b20d02c2
transformationMatrix(xhat, π/2, [1, 2, 3])

# ╔═╡ 631b3ec2-448b-11eb-2547-f93468ad8062
md"###### Adjoint Matrix"

# ╔═╡ e1e30074-4d54-11eb-3486-d352fd60e0f4
md"Computes the adjoint matrix of the transformation matrix `T`. Do not use a capital 'a'."

# ╔═╡ 4851772e-42aa-11eb-3259-ab30459bd137
md"###### Skew Matrix"

# ╔═╡ c90157b8-4d54-11eb-1c64-3ff71cc5fbd3
md"Computes the skew matrix of a vector or matrix `v` as defined in the book."

# ╔═╡ f3dc6268-42a8-11eb-0192-0b990bc0c84b
function skew(v)
	if length(v) == 6
		S = v
		ω = S[1:3]
		ν = S[4:6]
		return [skew(ω) ν; zeros(1, 3) 0]
	else
		ω = v
		return [0 -ω[3] ω[2]; ω[3] 0 -ω[1]; -ω[2] ω[1] 0]
	end
end

# ╔═╡ 63d3c80c-448b-11eb-0dc1-9b976df6e106
function adjoint(T::Array{Float64, 2})
	return [R(T) zeros(3, 3); skew(p(T))*R(T) R(T)]
end

# ╔═╡ 9af35988-448b-11eb-12f9-c7f3bee79403
adjoint(transformationMatrix(xhat, π/2, [1, 2, 3]))

# ╔═╡ d9c8c56e-42a9-11eb-3ac2-3fbe7f6991dd
skew([1, 2, 3, 4, 5, 6])

# ╔═╡ 353feda4-4d5a-11eb-27e8-11c2c259464f
md"#### Robots"

# ╔═╡ 3d9ecce0-4d5a-11eb-3da1-d9d48a6b15fa
struct Robot
	joints::Array{Array{Float64, 1}, 1}
	axis::Array{Array{Float64, 1}, 1}
	range::Float64
end

# ╔═╡ 993f82e4-4d58-11eb-03c2-0181e03d5121
md"3R robot as seen in the book p. 138."

# ╔═╡ a31ae364-4d5b-11eb-0c1b-4d52a89fd1e9
RRR = Robot([xhat, xhat, xhat], [yhat, yhat, yhat], 3)

# ╔═╡ 0bd4a532-4d59-11eb-0834-2578f2865d9a
md"KUKA LBR iiwa 7 d.o.f robot, as seen in task 4.1 and 5.1 in the exercise manual."

# ╔═╡ 1c554a76-4d5c-11eb-3532-17486d770041
begin
	KUKA = Robot([[0.0]], [[0.0]], 1)
	let L1 = 0.34, L2 = 0.4, L3 = 0.4, L4 = 0.15 
		KUKA = Robot([
			[0.0, 0.0, L1, 1.0],
			[0.0, 0.0, 0.0, 1.0],
			[0.0, 0.0, L2, 1.0],
			[0.0, 0.0, 0.0, 1.0],
			[0.0, 0.0, L3, 1.0],
			[0.0, 0.0, 0.0, 1.0],
			[0.0, 0.0, L4, 1.0]
		],
		[
			zhat,
			xhat,
			zhat,
			xhat,
			zhat,
			xhat,
			zhat
		], 1.5)
	end
end

# ╔═╡ 6545b464-4221-11eb-1dee-89c3d01e81de
md"#### Denavit Hartenberg parameters (Ch. 3.3)"

# ╔═╡ 6b1b4468-4d47-11eb-1492-e14006a0d20c
md"Returns the `transformation matrices` from {s} to {joint} for given joint `axis`, `joint` representations in their own reference frame and joint parameters `θ`. `
Tl` is the transformation matrix from {s} to {last joint} and `T` is the transformation matrix from {last joint} to {next joint}."

# ╔═╡ 5545204a-4d45-11eb-305e-0fb7d2133421
function transformationsDH(joints::Array{Array{Float64, 1}, 1}, axis::Array{Array{Float64, 1}, 1}, θ)
	
	ground::Array{Float64, 1} = [0.0, 0.0, 0.0, 1.0]
	T::Array{Float64, 2} = transformationMatrix(axis[1], θ[1], ground)
	transformations = [T]
	
	for index = 2:length(joints)
		Tl::Array{Float64, 2} = transformationMatrix(axis[index], θ[index], joints[index-1])
		T = T * Tl
		push!(transformations, T)
	end
	
	Tl = transformationMatrix(eye(3), joints[length(joints)])
	T = T * Tl
	push!(transformations, T)
	
	return transformations
end

# ╔═╡ 4ee62042-4d4c-11eb-0771-057f7ccfcaf2
transformationsDH(RRR.joints, RRR.axis, zeros(4))

# ╔═╡ 6e869bba-4221-11eb-191c-91a26a57e4e6
md"#### Plotting"

# ╔═╡ 631961c0-4d52-11eb-3200-27246f80c9c5
md"This function takes an array of transformation matrices `Ts` and plots a graph of the coordinate frames at each joint as well as the joints. This is set within the plot `limits`."

# ╔═╡ c10c9d9a-4d56-11eb-1aa3-0d8bc76ea222
md"This function takes arrays of `joints` represented in their reference frame, `axis` of rotation and parameters `θ` and plots a graph of the coordinate frames at each joint as well as the joints. This is set within the plot `limits`."

# ╔═╡ a3b36f62-4d56-11eb-39c9-6f20177d2cbe
plotRobot(R::Robot, θ) = plotRobot(transformationsDH(R.joints, R.axis, θ), (-R.range, R.range))

# ╔═╡ 42de89f8-4d54-11eb-0d8e-ad63a5b4ed9c
md"Displays the axis at a reference frame given by a transformation matrix `T`. Takes the parameter `plot` for the plot to draw in, this plot has to be actively displayed afterwards. Also takes the plot `limits` which dictate the length of the axis to make sure it is adequate."

# ╔═╡ 67f55590-3fa4-11eb-1904-5da71ba26e3a
function showAxis!(plot, T::Array{Float64, 2}, limits)
	l::Float64 = (limits[2] - limits[1]) / 20
	xAxis = T * [l, 0.0, 0.0, 1.0]
	yAxis = T * [0.0, l, 0.0, 1.0]
	zAxis = T * [0.0, 0.0, l, 1.0]
	origin = T * zero
	
	plot!(plot, [origin xAxis][1, :], [origin xAxis][2, :], [origin xAxis][3, :], color=:red)
	plot!(plot, [origin yAxis][1, :], [origin yAxis][2, :], [origin yAxis][3, :], color=:blue)
	plot!(plot, [origin zAxis][1, :], [origin zAxis][2, :], [origin zAxis][3, :], color=:green)
end

# ╔═╡ 8e4a9832-4d54-11eb-364b-f5df1ecebc55
md"Plots the joints of a robot which is represented by transformation matrices `Ts`. The joint coordinates are extracted from the matrices offset vector."

# ╔═╡ 135b0e0e-4d45-11eb-1870-87e23668d01b
function plotJoints(Ts::Array{Array{Float64, 2}, 1}, lims)
	l = length(Ts)
	xCord = [Ts[index][1, 4] for index = 1:l]
	yCord = [Ts[index][2, 4] for index = 1:l]
	zCord = [Ts[index][3, 4] for index = 1:l]
	plot(xCord, yCord, zCord, lims=lims, color=:black, linewidth=4, legend=false)
	scatter!(xCord, yCord, zCord, color=:white)
end

# ╔═╡ 7843cf14-4d4e-11eb-350d-0108281e1beb
function plotRobot(Ts::Array{Array{Float64, 2}, 1}, limits)
	l = length(Ts)
	p = plotJoints(Ts, limits)
	for index = 1:l
		showAxis!(p, Ts[index], limits)
	end
	p
end

# ╔═╡ 43f42af0-3f9e-11eb-202a-e7feee93095b
begin
	θ1Slider = @bind θ1 Slider(0:0.01:2π, show_value=true)
	θ2Slider = @bind θ2 Slider(0:0.01:2π, show_value=true)
	θ3Slider = @bind θ3 Slider(0:0.01:2π, show_value=true)
	θ4Slider = @bind θ4 Slider(0:0.01:2π, show_value=true)
	θ5Slider = @bind θ5 Slider(0:0.01:2π, show_value=true)
	θ6Slider = @bind θ6 Slider(0:0.01:2π, show_value=true)
	θ7Slider = @bind θ7 Slider(0:0.01:2π, show_value=true)
	
	md"""
	`θ1 = `$(θ1Slider)
	
	`θ2 = `$(θ2Slider)
	
	`θ3 = `$(θ3Slider)
	
	`θ4 = `$(θ4Slider)
	
	`θ5 = `$(θ5Slider)
	
	`θ6 = `$(θ6Slider)
	
	`θ7 = `$(θ7Slider)
	"""
end

# ╔═╡ 540b2a76-3f86-11eb-1426-c353435545d2
plotRobot(RRR, [θ1, θ2, θ3])

# ╔═╡ 0b02c060-4d5d-11eb-0b56-87148bf115ec
plotRobot(KUKA, [θ1, θ2, θ3, θ4, θ5, θ6, θ7])

# ╔═╡ 8c201bbe-4489-11eb-3960-e34b74493aad
md"#### Manipulator Jacobian (Ch. 5)"

# ╔═╡ d2dcea70-4489-11eb-0aca-57d973131625
begin
	S::Array{Float64, 2} = [0.0 0.0; 0.0 0.0]
	B::Array{Float64, 2} = [0.0 0.0; 0.0 0.0]
	let L1 = 0.34, L2 = 0.4, L3 = 0.4, L4 = 0.15
		S = Transpose([
			0 0 1 0 0 0;
			1 0 0 0 L1 0;
			0 0 1 0 0 0;
			1 0 0 0 L1+L2 0;
			0 0 1 0 0 0;
			1 0 0 0 L1+L2+L3 0;
			0 0 1 0 0 0
		])

		B = Transpose([
			0 0 1 0 0 0;
			1 0 0 0 -L4-L3-L2 0;
			0 0 1 0 0 0;
			1 0 0 0 -L4-L3 0;
			0 0 1 0 0 0;
			1 0 0 0 -L4 0;
			0 0 1 0 0 0
		])
	end
end

# ╔═╡ acd33468-4489-11eb-1d06-17d278e4b389
md"###### Space Jacobian from screws"

# ╔═╡ 5f1f14da-4d53-11eb-25d6-3b5fcf15b3aa
md"Computes the Jacobian in {s} using the `screws` and `θ`."

# ╔═╡ 95052d82-4489-11eb-0406-a9efde904d63
function spaceJacobian(screws::Array{Float64, 2}, θ::Array{Float64, 1})
	n = length(θ)
	Jacobian::Array{Float64, 2} = eye(6, n)
	for index in 1:n
		Jacobian[:, index] = spaceJacobian(screws, θ, index)
	end
	return Jacobian
end

# ╔═╡ 468158b6-4d53-11eb-155d-f9ba5c3efbfb
md"Returns the `i`th column of the jacobian."

# ╔═╡ abd4cf40-4d30-11eb-183a-835ccc945984
function spaceJacobian(screws::Array{Float64, 2}, θ::Array{Float64, 1}, i::Int64)
	Jsn = eye(4)
	for index = 1:i-1
		Jsn = Jsn * exp(skew(screws[:, index]) * θ[index])
	end
	return adjoint(Jsn) * screws[:, i]
end

# ╔═╡ 5c346f58-448a-11eb-2ca0-ab84dde13f94
spaceJacobian(S, zeros(7))

# ╔═╡ bc097d8c-4d41-11eb-05d3-a9c021debf29
md"###### Body Jacobian from screws"

# ╔═╡ 0fb20fb2-4d60-11eb-298f-6d1a849b8d79
md"Computes the Jacobian in {b} using the `screws` and `θ`."

# ╔═╡ cf0b6814-4d41-11eb-3f3b-5fd582554343
function bodyJacobian(screws::Array{Float64, 2}, θ::Array{Float64, 1})
	n = length(θ)
	Jacobian::Array{Float64, 2} = eye(6, n)
	for index in 1:n
		Jacobian[:, index] = bodyJacobian(screws, θ, index)
	end
	return Jacobian
end

# ╔═╡ 172cc686-4d60-11eb-2799-ad1fc4a2712c
md"Returns the `i`th column of the jacobian."

# ╔═╡ e5e0f806-4d41-11eb-3d5a-075ebee9ab8b
function bodyJacobian(screws::Array{Float64, 2}, θ::Array{Float64, 1}, i::Int64)
	Jsn = eye(4)
	for index = length(θ):-1:i+1
		Jsn = Jsn * exp(-skew(screws[:, index]) * θ[index])
	end
	return adjoint(Jsn) * screws[:, i]
end

# ╔═╡ f15b295a-4d5e-11eb-24a0-7944dd834bd6
bodyJacobian(B, zeros(7))

# ╔═╡ Cell order:
# ╟─54276a22-3f80-11eb-2538-e18635928ecf
# ╠═36844d82-3f80-11eb-244a-83efb39d115d
# ╟─350a9f42-3f85-11eb-06ae-2556cc9aa4fd
# ╟─fb3fedce-4d55-11eb-2879-b755a16cde4e
# ╠═327e9936-3f85-11eb-1bf9-d1a2def49b41
# ╟─c7f5bac4-4da1-11eb-333f-63859c65f8f7
# ╟─83fd510c-42ab-11eb-347c-4f9625869adc
# ╟─dcdfa036-4d55-11eb-1026-1763c6107e31
# ╠═8d4bdea4-42ab-11eb-1f3a-e5cea01bcf3f
# ╟─f2a1fbbc-4d55-11eb-0ff1-41e28c57ad71
# ╠═faa7d0d0-4d36-11eb-38e8-371dc6654639
# ╟─f0952168-3f82-11eb-37f2-b17628f1044c
# ╟─9be3e8c6-4d55-11eb-1949-43d88b18b0db
# ╠═cd374f48-3f82-11eb-26f6-8dfbc7216940
# ╟─6e303e66-4d55-11eb-3c6f-59cfcbf9211c
# ╠═5e89ba96-3f81-11eb-0a35-e9d3928427b8
# ╠═30385498-3f83-11eb-1584-792e07e24dc7
# ╟─0829b69a-3f83-11eb-3f65-0f5ab0b35aa9
# ╟─176e3a1a-4d55-11eb-0dd2-27d5cdac3b20
# ╠═152049ea-3f83-11eb-18dd-8d667de8625d
# ╟─3abe0ebe-4d55-11eb-05f8-35ce4ffefd66
# ╠═09009752-3f86-11eb-1e79-3f70a93d7575
# ╟─c2826700-4d55-11eb-19cf-7b814acfca64
# ╠═0d673ae0-42ac-11eb-05ec-7798c0839217
# ╟─cf770c18-4d55-11eb-362f-a138a8b6c353
# ╠═1db602d0-42ac-11eb-0992-735dc410de36
# ╠═dee42ada-3f83-11eb-2c39-83a2b20d02c2
# ╟─631b3ec2-448b-11eb-2547-f93468ad8062
# ╟─e1e30074-4d54-11eb-3486-d352fd60e0f4
# ╠═63d3c80c-448b-11eb-0dc1-9b976df6e106
# ╠═9af35988-448b-11eb-12f9-c7f3bee79403
# ╟─4851772e-42aa-11eb-3259-ab30459bd137
# ╟─c90157b8-4d54-11eb-1c64-3ff71cc5fbd3
# ╠═f3dc6268-42a8-11eb-0192-0b990bc0c84b
# ╠═d9c8c56e-42a9-11eb-3ac2-3fbe7f6991dd
# ╟─353feda4-4d5a-11eb-27e8-11c2c259464f
# ╠═3d9ecce0-4d5a-11eb-3da1-d9d48a6b15fa
# ╟─993f82e4-4d58-11eb-03c2-0181e03d5121
# ╟─a31ae364-4d5b-11eb-0c1b-4d52a89fd1e9
# ╟─0bd4a532-4d59-11eb-0834-2578f2865d9a
# ╟─1c554a76-4d5c-11eb-3532-17486d770041
# ╟─6545b464-4221-11eb-1dee-89c3d01e81de
# ╟─6b1b4468-4d47-11eb-1492-e14006a0d20c
# ╠═5545204a-4d45-11eb-305e-0fb7d2133421
# ╠═4ee62042-4d4c-11eb-0771-057f7ccfcaf2
# ╟─6e869bba-4221-11eb-191c-91a26a57e4e6
# ╟─631961c0-4d52-11eb-3200-27246f80c9c5
# ╠═7843cf14-4d4e-11eb-350d-0108281e1beb
# ╟─c10c9d9a-4d56-11eb-1aa3-0d8bc76ea222
# ╠═a3b36f62-4d56-11eb-39c9-6f20177d2cbe
# ╟─42de89f8-4d54-11eb-0d8e-ad63a5b4ed9c
# ╠═67f55590-3fa4-11eb-1904-5da71ba26e3a
# ╟─8e4a9832-4d54-11eb-364b-f5df1ecebc55
# ╠═135b0e0e-4d45-11eb-1870-87e23668d01b
# ╠═540b2a76-3f86-11eb-1426-c353435545d2
# ╟─43f42af0-3f9e-11eb-202a-e7feee93095b
# ╠═0b02c060-4d5d-11eb-0b56-87148bf115ec
# ╟─8c201bbe-4489-11eb-3960-e34b74493aad
# ╠═d2dcea70-4489-11eb-0aca-57d973131625
# ╟─acd33468-4489-11eb-1d06-17d278e4b389
# ╟─5f1f14da-4d53-11eb-25d6-3b5fcf15b3aa
# ╠═95052d82-4489-11eb-0406-a9efde904d63
# ╟─468158b6-4d53-11eb-155d-f9ba5c3efbfb
# ╠═abd4cf40-4d30-11eb-183a-835ccc945984
# ╠═5c346f58-448a-11eb-2ca0-ab84dde13f94
# ╟─bc097d8c-4d41-11eb-05d3-a9c021debf29
# ╟─0fb20fb2-4d60-11eb-298f-6d1a849b8d79
# ╠═cf0b6814-4d41-11eb-3f3b-5fd582554343
# ╟─172cc686-4d60-11eb-2799-ad1fc4a2712c
# ╠═e5e0f806-4d41-11eb-3d5a-075ebee9ab8b
# ╠═f15b295a-4d5e-11eb-24a0-7944dd834bd6
