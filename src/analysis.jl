import LazySets: convex_hull

function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	tracking_data = npzread(filename)[1, :, 1:TRACKING_DIM, :]
	#TODO: warn when contains NaNs?
	# Get the number of robots and the number of timesteps
	n_robots = size(tracking_data, ROBOTS)

	# Add heading vector from orientation to data
	heading_vector_xs = cos.(tracking_data[:, θ:θ, :])
	heading_vector_ys = sin.(tracking_data[:, θ:θ, :])
	velocity_xs = cat(zeros(n_robots, 1, 1), diff(tracking_data[:, X:X, :]; dims=T); dims=T)
	velocity_ys = cat(zeros(n_robots, 1, 1), diff(tracking_data[:, Y:Y, :]; dims=T); dims=T)
	velocity_magnitude = sqrt.(velocity_xs .^ 2 .+ velocity_ys .^ 2)
	acceleration_xs = cat(diff(velocity_xs; dims=T), zeros(n_robots, 1, 1); dims=T)
	acceleration_ys = cat(diff(velocity_ys; dims=T), zeros(n_robots, 1, 1); dims=T)
	acceleration_magnitude = sqrt.(acceleration_xs .^ 2 .+ acceleration_ys .^ 2)
	angular_velocity = cat(
		zeros(n_robots, 1, 1),
		[
			abs(d) < 2pi - abs(d) ? d : -sign(d) * (2pi - abs(d)) for
			d in diff(tracking_data[:, θ:θ, :]; dims=T)
		];
		dims=T,
	)
	angular_acceleration = cat(
		diff(angular_velocity; dims=T), zeros(n_robots, 1, 1); dims=T
	)

	tracking_data = cat(
		tracking_data,
		heading_vector_xs,
		heading_vector_ys,
		velocity_xs,
		velocity_ys,
		velocity_magnitude,
		acceleration_xs,
		acceleration_ys,
		acceleration_magnitude,
		angular_velocity,
		angular_acceleration;
		dims=PROPERTIES,
	)

	# precalculate all metrics for all timesteps (currently no clustering)
	#TODO: more outsourcing to metrics.jl?
	polarisation = swarm_polarisation.(eachslice(tracking_data; dims=T))
	rotational_order = swarm_rotational_order.(eachslice(tracking_data; dims=T))
	distmats = [
		pairwise(Euclidean(), permutedims(s)) for
		s in eachslice(tracking_data[:, [X, Y], :]; dims=T)
	]
	mean_interindividual_distance = mean.(distmats)
	surrounding_polygon =
		convex_hull.(
			collect(eachslice(s; dims=ROBOTS)) for
			s in eachslice(tracking_data[:, [X, Y], :]; dims=T)
		)
	center_of_mass = [
		mean(eachslice(s; dims=ROBOTS)) for
		s in eachslice(tracking_data[:, [X, Y], :]; dims=T)
	]
	findmaxdist = [findmax(m) for m in distmats]
	diameter = [f[1] for f in findmaxdist]
	furthest = [f[2] for f in findmaxdist]
	maxmindist = [
		maximum(minimum.(eachcol(m + diagm(Inf * ones(size(m, 1)))))) for m in distmats
	]
	minmaxdist = [minimum(maximum.(eachcol(m))) for m in distmats]
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

	metrics = Dict(
		"Polarisation" => polarisation,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
		"Max Min IID" => maxmindist,
		"Min Max IID" => minmaxdist,
	)
	geometry = Dict( #TODO: rename
		"Surrounding Polygon" => surrounding_polygon,
		"Center of Mass" => center_of_mass,
		"Furthest Robots" => furthest,
		"Distance Matrices" => distmats,
	)

	return SwarmData(
		tracking_data[:, 1:TRACKING_DIM, :],
		tracking_data[:, TRACKING_DIM:end, :],
		metrics,
		geometry,
	)
end

function analyse_wall(filename)
	wd = npzread(filename)[1, 1, [X, Y], :]
	# Remove columns with NaN values from wall_data
	wd = wd[:, .!isnan.(wd[1, :])]
	wd = wd[:, .!isnan.(wd[2, :])]
	return wd
end

function process_collisions(filename)
	collision_dict = JSON.parse(String(read(filename)))["0"]
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
