import LazySets: convex_hull
using Clustering
using DSP
using StatsBase

# choose weighting of the total variation programmaticially like in
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7899139/
# but with the power-weighted mean frequency instead of choosing it manually
function fitted_gamma(data, fs)
	periodograms = periodogram.(eachslice(data; dims=(ROBOTS, PROPERTIES)), fs=fs)
	mean_power = mean(vec([pg.power for pg in periodograms]))
	mean_frequency = mean(periodograms[1].freq, weights(mean_power))
	return exp(-1.6 * log(mean_frequency) + 0.71 * log(fs) - 5.1)
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
	fs = 1 ÷ mean(diff(robot_data[:, T:T, :]; dims=TIME))
	# TODO: tvr diff
	velocity_xs, velocity_zs =
		[diff(robot_data[:, i:i, :]; dims=TIME) for i in (X, Z)] .* fs
	velocity_magnitude = sqrt.(velocity_xs .^ 2 .+ velocity_zs .^ 2)
	acceleration_xs, acceleration_zs =
		[diff(vs; dims=TIME) for vs in (velocity_xs, velocity_zs)] .* fs
	acceleration_magnitude = sqrt.(acceleration_xs .^ 2 .+ acceleration_zs .^ 2)

	robot_data = cat(
		robot_data,
		heading_vector_xs,
		heading_vector_zs,
		velocity_xs,
		velocity_zs,
		velocity_magnitude,
		acceleration_xs,
		acceleration_zs,
		acceleration_magnitude;
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
	acceleration_correlation = [
		mean(dot(a, s[j, :]) for (i, a) in enumerate(eachrow(s)) for j in 1:(i - 1)) for
		s in eachslice(robot_data[:, [ACCX, ACCZ], :]; dims=TIME)
	]
	cosine_dist =
		(
			[
				pairwise(Cosine(), permutedims(s)) for
				s in eachslice(robot_data[:, [X, Z], :]; dims=TIME)
			] .+ 1
		) ./ 2
	velocity_dist = [
		[min(v[i] / v[j], v[j] / v[i]) for i in eachindex(v), j in 1:(i - 1)] for
		v in eachslice(velocity_magnitude; dims=TIME)
	]
	min_dist = minimum(stack(distmats))
	visited_diameter = sqrt(
		sum([reduce(-, extrema(robot_data[:, i, :]))^2 for i in (X, Z)])
	)
	clusterings =
		hclust.(
			(distmats .- min_dist) ./ (visited_diameter - min_dist) .*
			(visited_diameter / min_dist) .+ cosine_dist .+ velocity_dist, #TODO add weights
			branchorder=:barjoseph,
		)
	single_cluster_thresholds = [maximum(c.heights) for c in clusterings]
	two_cluster_thresholds = [sort(c.heights)[end - 1] for c in clusterings]

	metrics = Dict(
		"Polarisation" => polarisation,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
		"Max Min IID" => maxmindist,
		"Correlation of Acceleration" => acceleration_correlation,
		"1 to 2 Clusters Threshold" => single_cluster_thresholds,
		"2 to 3 Clusters Threshold" => two_cluster_thresholds,
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
