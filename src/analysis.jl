import LazySets: convex_hull
using Clustering
using StatsBase

function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	robot_data = npzread(filename)[1, :, 1:TRACKING_DIM, :]
	#TODO: warn when contains NaNs?

	# Add heading vector from orientation to data
	robot_data[:, θ, :] = mod2pi.(robot_data[:, θ, :])
	heading_vector_xs = cos.(robot_data[:, θ:θ, :])
	heading_vector_zs = sin.(robot_data[:, θ:θ, :])

	robot_data = cat(robot_data, heading_vector_xs, heading_vector_zs; dims=PROPERTIES)

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

	cosine_distance_matrices =
		[
			pairwise(CosineDist(), permutedims(s)) for
			s in eachslice(robot_data[:, [HVX, HVZ], :]; dims=TIME)
		] ./ 2
	visited_diameter = sqrt(
		sum([reduce(-, extrema(robot_data[:, i, :]))^2 for i in (X, Z)])
	)
	dissimilarity_matrices = [
		1 .- ((1 .- dm ./ visited_diameter) .* (1 .- cdm)) .^ (1 / 3) for
		(dm, cdm) in zip(distance_matrices, cosine_distance_matrices)
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
