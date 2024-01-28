import LazySets: convex_hull
using Clustering

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
	velocity_xs = cat(zeros(n_robots, 1, 1), diff(robot_data[:, X:X, :]; dims=T); dims=T)
	velocity_zs = cat(zeros(n_robots, 1, 1), diff(robot_data[:, Z:Z, :]; dims=T); dims=T)
	velocity_magnitude = sqrt.(velocity_xs .^ 2 .+ velocity_zs .^ 2)
	acceleration_xs = cat(diff(velocity_xs; dims=T), zeros(n_robots, 1, 1); dims=T)
	acceleration_zs = cat(diff(velocity_zs; dims=T), zeros(n_robots, 1, 1); dims=T)
	acceleration_magnitude = sqrt.(acceleration_xs .^ 2 .+ acceleration_zs .^ 2)
	angular_velocity = cat( #TODO: broken until input is fixed
		zeros(n_robots, 1, 1),
		[ #TODO reshape?
			abs(d) < 2pi - abs(d) ? d : -sign(d) * (2pi - abs(d)) for
			d in diff(robot_data[:, θ:θ, :]; dims=T)
		];
		dims=T,
	)
	angular_acceleration = cat(
		diff(angular_velocity; dims=T), zeros(n_robots, 1, 1); dims=T
	)

	robot_data = cat(
		robot_data,
		heading_vector_xs,
		heading_vector_zs,
		velocity_xs,
		velocity_zs,
		velocity_magnitude,
		acceleration_xs,
		acceleration_zs,
		acceleration_magnitude,
		angular_velocity,
		angular_acceleration;
		dims=PROPERTIES,
	)

	# precalculate all metrics for all timesteps (currently no clustering)
	#TODO: more outsourcing to metrics.jl?
	polarisation = swarm_polarisation.(eachslice(robot_data; dims=T))
	rotational_order = swarm_rotational_order.(eachslice(robot_data; dims=T))
	distmats = [
		pairwise(Euclidean(), permutedims(s)) for
		s in eachslice(robot_data[:, [X, Z], :]; dims=T)
	]
	mean_interindividual_distance = mean.(distmats) .* (1 .+ 1 ./ (size.(distmats, 1) .- 1))
	surrounding_polygon =
		convex_hull.(
			collect(eachslice(s; dims=ROBOTS)) for
			s in eachslice(robot_data[:, [X, Z], :]; dims=T)
		)
	center_of_mass = [
		mean(eachslice(s; dims=ROBOTS)) for s in eachslice(robot_data[:, [X, Z], :]; dims=T)
	]
	findmaxdist = [findmax(m) for m in distmats]
	diameter = [f[1] for f in findmaxdist]
	furthest = [[f[2][1], f[2][2]] for f in findmaxdist]
	maxmindist = [
		maximum(minimum.(eachcol(m + diagm(Inf * ones(size(m, 1)))))) for m in distmats
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
	# acceleration_correlation = [
	# 	mean(dot(a, s[j, :]) for (i, a) in enumerate(eachrow(s)) for j in 1:(i - 1)) for
	# 	s in eachslice(smoothed_acceleration; dims=T)
	# ]

	clusterings = hclust.(distmats, branchorder=:barjoseph)
	single_cluster_thresholds = [maximum(c.heights) for c in clusterings]

	metrics = Dict(
		"Polarisation" => polarisation,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
		"Max Min IID" => maxmindist,
		"Single Cluster Threshold" => single_cluster_thresholds,
	)
	derived = Dict( #TODO: rename
		"Convex Hull" => surrounding_polygon,
		"Center of Mass" => center_of_mass,
		"Furthest Robots" => furthest,
		"Distance Matrices" => distmats,
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
