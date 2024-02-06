import LazySets: convex_hull
using Clustering
using DSP
using StatsBase

function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	robot_data = npzread(filename)[1, :, 1:TRACKING_DIM, :]
	#TODO: warn when contains NaNs?
	# Get the number of robots and the number of timesteps
	n_robots = size(robot_data, ROBOTS)

	# Add heading vector from orientation to data
	robot_data[:, θ, :] = mod2pi.(robot_data[:, θ, :])
	heading_vector_xs = cos.(robot_data[:, θ:θ, :])
	heading_vector_zs = sin.(robot_data[:, θ:θ, :])
    #TODO remove after final rejection
    # choose weighting of the total variation programmaticially like in
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7899139/
	# fs = round(Int, 1 / mean(diff(robot_data[:, T:T, :]; dims=TIME)))
	# velocity_xs, velocity_zs =
	# 	[
	# 		cat(zeros(n_robots), diff(robot_data[:, i:i, :]; dims=TIME); dims=TIME) for
	# 		i in (X, Z)
	# 	] .* fs
	# velocity_xs, velocity_zs = [
	# 	permutedims(
	# 		stack([
	# 			filt(
	# 				digitalfilter(Lowpass(0.25; fs=fs), Butterworth(1)),
	# 				tvdiff(
	# 					s,
	# 					40,
	# 					exp(-1.6 * log(0.25) + 0.71 * log(30) - 5.1);
	# 					scale="large",
	# 					dx=1 / 30,
	# 				),
	# 			) for s in eachslice(robot_data[:, i:i, :]; dims=(ROBOTS, PROPERTIES))
	# 		]),
	# 		(2, 3, 1),
	# 	) for i in (X, Z)
	# ]
	# velocity_magnitude = sqrt.(velocity_xs .^ 2 .+ velocity_zs .^ 2)
	# acceleration_xs, acceleration_zs =
	# 	[
	# 		cat(diff(vs; dims=TIME), zeros(n_robots); dims=TIME) for
	# 		vs in (velocity_xs, velocity_zs)
	# 	] .* fs
	# acceleration_xs, acceleration_zs = [
	# 	permutedims(
	# 		stack([
	# 			filt(
	# 				digitalfilter(Lowpass(0.25; fs=fs), Butterworth(1)),
	# 				tvdiff(
	# 					s,
	# 					40,
	# 					exp(-1.6 * log(0.25) + 0.71 * log(30) - 5.1);
	# 					scale="large",
	# 					dx=1 / 30,
	# 				),
	# 			) for s in eachslice(vs; dims=(ROBOTS, PROPERTIES))
	# 		]),
	# 		(2, 3, 1),
	# 	) for vs in (velocity_xs, velocity_zs)
	# ]

	# acceleration_magnitude = sqrt.(acceleration_xs .^ 2 .+ acceleration_zs .^ 2)

	robot_data = cat(
		robot_data,
		heading_vector_xs,
		heading_vector_zs,
		# velocity_xs, #TODO remove after final rejection
		# velocity_zs,
		# velocity_magnitude,
		# acceleration_xs,
		# acceleration_zs,
		# acceleration_magnitude;
		dims=PROPERTIES,
	)

	# precalculate all metrics for all timesteps (currently no clustering)
	#TODO: more outsourcing to metrics.jl?
	polarisation = swarm_polarisation.(eachslice(robot_data; dims=TIME))
	rotational_order = swarm_rotational_order.(eachslice(robot_data; dims=TIME))
	distance_matrices = [
		pairwise(Euclidean(), permutedims(s)) for
		s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
	]
	mean_interindividual_distance =
		mean.(distance_matrices) .* (1 .+ 1 ./ (size.(distance_matrices, 1) .- 1))
	surrounding_polygon =
		convex_hull.(
			collect(eachslice(s; dims=ROBOTS)) for
			s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
		)
	center_of_mass = [
		mean(eachslice(s; dims=ROBOTS)) for
		s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
	]
	findmaxdist = [findmax(m) for m in distance_matrices]
	diameter = [f[1] for f in findmaxdist]
	furthest = [[f[2][1], f[2][2]] for f in findmaxdist]
	maxmindist = [
		maximum(minimum.(eachcol(m + diagm(Inf * ones(size(m, 1)))))) for
		m in distance_matrices
	]
	area = [ #TODO: outsource
		0.5 * abs( # shoelace formula for area of polygon (https://en.wikipedia.org/wiki/Shoelace_formula)
			sum(
				p[i][1] * p[i % length(p) + 1][2] - p[i % length(p) + 1][1] * p[i][2] for
				i in eachindex(p)
			),
		) for p in surrounding_polygon
	]
	roundness = [
		4pi * A / sum(colwise(Euclidean(), stack(p), stack(vcat(p[2:end], [p[1]]))))^2 for
		(p, A) in zip(surrounding_polygon, area)
	]
    #TODO remove after final rejection
	# acceleration_moving_average = map(x -> conv(ones(61) / 61, x), acceleration_magnitude)
	# acceleration_correlation = [
	# 	mean(dot(a, s[j, :]) for (i, a) in enumerate(eachrow(s)) for j in 1:(i - 1)) for
	# 	s in eachslice(robot_data[:, [ACCX, ACCZ], :]; dims=TIME)
	# ]
	cosine_distance_matrices =
		[
			pairwise(CosineDist(), permutedims(s)) for
			s in eachslice(robot_data[:, [HVX, HVZ], :]; dims=TIME)
		] ./ 2
	# velocity_similarity_matrices = [
	# 	[(min(u, w) + eps(Float64)) / (max(u, w) + eps(Float64)) for u in v, w in v] for
	# 	v in eachslice(velocity_magnitude; dims=(TIME, PROPERTIES))
	# ]
	visited_diameter = sqrt(
		sum([reduce(-, extrema(robot_data[:, i, :]))^2 for i in (X, Z)])
	)
	dissimilarity_matrices = [
		1 .- ((1 .- dm ./ visited_diameter) .* (1 .- cdm)) .^ (1 / 3) for (dm, cdm, vsm) in
		zip(distance_matrices, cosine_distance_matrices, velocity_similarity_matrices)
	]
	clusterings = hclust.(dissimilarity_matrices, branchorder=:barjoseph)
	single_cluster_thresholds = [log(maximum(c.heights)) for c in clusterings]

	metrics = Dict(
		"Polarisation" => polarisation,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
		"Max Min IID" => maxmindist,
		"Log Last Merge Threshold" => single_cluster_thresholds,
		# "Corr. of Acceleration" => acceleration_correlation,
		# "VelX" => velocity_xs[1, 1, :],
		# "AccX" => acceleration_xs[1, 1, :],
		# "Vel" => velocity_magnitude[1, 1, :],
		# "Acc" => acceleration_magnitude[1, 1, :],
	)
	derived = Dict( #TODO: rename
		"Convex Hull" => surrounding_polygon,
		"Center of Mass" => center_of_mass,
		"Furthest Robots" => furthest,
		"Distance Matrices" => distance_matrices,
	)

	return SwarmData(robot_data, metrics, derived, clusterings)
end

function analyse_wall(filename)
	wd = npzread(filename)[1, 1, [X, Z], :]
	# Remove columns with NaN values from wall_data
	wd = wd[:, .!isnan.(wd[1, :])]
	wd = wd[:, .!isnan.(wd[2, :])]
	return wd
end

function process_collisions(filename)
	collision_dict = JSON3.read(filename, Dict)["0"]
	collisions = falses(
		length(keys(collision_dict)), maximum(vcat(values(collision_dict)...))
	)
	for (robot, collision_list) in collision_dict
		robot_row = parse(Int, robot) + 1
		for collision in collision_list
			collisions[robot_row, collision] = true
		end
	end
	return collisions
end
