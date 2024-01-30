import LazySets: convex_hull
using Clustering
using DSP

function apply_lowpass(data, cutoff, order, dts)
	return permutedims(
		filt(
			digitalfilter(Lowpass(cutoff; fs=1 ÷ mean(dts)), Butterworth(order)),
			permutedims(data .- data[:, :, 1], (3, 1, 2)),
		),
		(2, 3, 1),
	) .+ data[:, :, 1]
end

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
	dts =
		mean.(eachslice(diff(robot_data[:, T:T, :]; dims=TIME); dims=(ROBOTS, PROPERTIES)))
	xs_lowpass, zs_lowpass = [
		apply_lowpass(robot_data[:, i:i, :], 0.16, 1, dts) for i in (X, Z)
	]

	velocity_xs, velocity_zs = [
		cat(
			zeros(n_robots, 1, 1),
			apply_lowpass(diff(pos_lp; dims=TIME) ./ dts, 0.16, 1, dts);
			dims=TIME,
		) for pos_lp in (xs_lowpass, zs_lowpass)
	]
	velocity_magnitude = sqrt.(velocity_xs .^ 2 .+ velocity_zs .^ 2)
	acceleration_xs, acceleration_zs = [
		cat(
			apply_lowpass(diff(vel_lp; dims=TIME) ./ dts, 0.16, 1, dts),
			zeros(n_robots, 1, 1);
			dims=TIME,
		) for vel_lp in (velocity_xs, velocity_zs)
	]
	acceleration_magnitude = sqrt.(acceleration_xs .^ 2 .+ acceleration_zs .^ 2)
	angular_velocity = cat(
		zeros(n_robots, 1, 1),
		map(
			d -> abs(d) < pi ? d : -sign(d) * (2pi - abs(d)),
			diff(robot_data[:, θ:θ, :]; dims=TIME),
		) ./ dts;
		dims=TIME,
	)
	angular_velocity_lowpass = permutedims(
		filt(
			digitalfilter(Lowpass(0.05; fs=30), Butterworth(1)),
			permutedims(angular_velocity[1:10, :, :], (3, 1, 2)),
		),
		(2, 3, 1),
	)
	angular_acceleration = cat(
		diff(angular_velocity_lowpass; dims=TIME), zeros(n_robots, 1, 1); dims=TIME
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
	polarisation = swarm_polarisation.(eachslice(robot_data; dims=TIME))
	rotational_order = swarm_rotational_order.(eachslice(robot_data; dims=TIME))
	distmats = [
		pairwise(Euclidean(), permutedims(s)) for
		s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
	]
	mean_interindividual_distance = mean.(distmats) .* (1 .+ 1 ./ (size.(distmats, 1) .- 1))
	surrounding_polygon =
		convex_hull.(
			collect(eachslice(s; dims=ROBOTS)) for
			s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
		)
	center_of_mass = [
		mean(eachslice(s; dims=ROBOTS)) for
		s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
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
	# 	s in eachslice(robot_data[:, [ACCX, ACCZ], :]; dims=TIME)
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
		"Correlation of Acceleration" => acceleration_correlation,
		"Single Cluster Threshold" => single_cluster_thresholds,
		"X_lowpass" => xs_lowpass[1, 1, :],
		"VEL" => velocity_xs[1, 1, :],
		"Acc" => acceleration_xs[1, 1, :],
		"Angle" => robot_data[1, θ, :],
		"dAngle" => angular_velocity[1, 1, :],
		"dAngle_lowpass" => angular_velocity_lowpass[1, 1, :],
		"ddAngle" => angular_acceleration[1, 1, :],
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
