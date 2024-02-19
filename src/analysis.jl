import LazySets: convex_hull
using Clustering
using StatsBase

function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	agent_data = npzread(filename)[1, :, 1:TRACKING_DIM, :]

	# Add heading vector from orientation to data
	agent_data[:, θ, :] = mod2pi.(agent_data[:, θ, :])
	heading_vector_xs = cos.(agent_data[:, θ:θ, :])
	heading_vector_zs = sin.(agent_data[:, θ:θ, :])

	agent_data = cat(agent_data, heading_vector_xs, heading_vector_zs; dims=PROPERTIES)

	# precalculate all metrics for all timesteps (currently no clustering)
	center_of_mass = [
		mean(eachslice(s; dims=AGENTS)) for
		s in eachslice(agent_data[:, [X, Z], :]; dims=TIME)
	]
	polarization = swarm_polarization.(eachslice(agent_data; dims=TIME))
	rotational_order =
		swarm_rotational_order.(eachslice(agent_data; dims=TIME), center_of_mass)
	distance_matrices = [
		pairwise(Euclidean(), permutedims(s)) for
		s in eachslice(agent_data[:, [X, Z], :]; dims=TIME)
	]
	mean_interindividual_distance =
		mean.(distance_matrices) .* (1 .+ 1 ./ (size.(distance_matrices, 1) .- 1))
	convex_hulls =
		convex_hull.(
			collect(eachslice(s; dims=AGENTS)) for
			s in eachslice(agent_data[:, [X, Z], :]; dims=TIME)
		)

	findmaxdist = [findmax(m) for m in distance_matrices]
	diameter = [f[1] for f in findmaxdist]
	furthest = [[f[2][1], f[2][2]] for f in findmaxdist]
	maxmindist = [
		maximum(minimum.(eachcol(m + diagm(Inf * ones(size(m, 1)))))) for
		m in distance_matrices
	]
	area = shoelace_area.(convex_hulls)
	roundness = [
		4pi * A / sum(colwise(Euclidean(), stack(p), stack(circshift(p, -1))))^2 for
		(p, A) in zip(convex_hulls, area)
	]
	mean_mindist = [
		mean(minimum.(eachcol(m + diagm(Inf * ones(size(m, 1)))))) for
		m in distance_matrices
	]

	cosine_distance_matrices =
		[
			pairwise(CosineDist(), permutedims(s)) for
			s in eachslice(agent_data[:, [HVX, HVZ], :]; dims=TIME)
		] ./ 2
	visited_diameter = sqrt(
		sum([reduce(-, extrema(agent_data[:, i, :]))^2 for i in (X, Z)])
	)
	dissimilarity_matrices = [
		1 .- sqrt.((1 .- dm ./ visited_diameter) .* (1 .- cdm)) for
		(dm, cdm) in zip(distance_matrices, cosine_distance_matrices)
	]
	clusterings = hclust.(dissimilarity_matrices, branchorder=:barjoseph)
	single_cluster_thresholds = [log(maximum(c.heights)) for c in clusterings]

	metrics = Dict(
		"Polarization" => polarization,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
		"Max Min IID" => maxmindist,
        "Mean Min IID" => mean_mindist,
		"Log Last Merge Threshold" => single_cluster_thresholds,
	)
	derived = Dict(
		"Convex Hull" => convex_hulls,
		"Center of Mass" => center_of_mass,
		"Furthest Robots" => furthest,
		"Distance Matrices" => distance_matrices,
	)

	return SwarmData(agent_data, metrics, derived, clusterings)
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
	for (agent, collision_list) in collision_dict
		agent_row = parse(Int, agent) + 1
		for collision in collision_list
			collisions[agent_row, collision] = true
		end
	end
	return collisions
end
