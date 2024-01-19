import LazySets: convex_hull

function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	tracking_data = npzread(filename)[1, :, 1:TRACKING_DIM, :]
	#TODO: warn when contains NaNs?
	# Get the number of robots and the number of timesteps
	n_robots = size(tracking_data, ROBOTS)
	n_timesteps = size(tracking_data, T)

	# Add heading vector from orientation to data
	heading_vector_xs = reshape(cos.(tracking_data[:, θ, :]), n_robots, 1, n_timesteps)
	heading_vector_ys = reshape(sin.(tracking_data[:, θ, :]), n_robots, 1, n_timesteps)
	velocity_xs = reshape(
		cat(
			zeros(n_robots), diff(tracking_data[:, X, :]; dims=PROPERTIES); dims=PROPERTIES
		),
		n_robots,
		1,
		n_timesteps,
	)
	velocity_ys = reshape(
		cat(
			zeros(n_robots), diff(tracking_data[:, Y, :]; dims=PROPERTIES); dims=PROPERTIES
		),
		n_robots,
		1,
		n_timesteps,
	)
	velocity_magnitude = reshape(
		sqrt.(velocity_xs .^ 2 .+ velocity_ys .^ 2), n_robots, 1, n_timesteps
	)
	acceleration_xs = reshape(
		cat(diff(velocity_xs; dims=T), zeros(n_robots); dims=T), n_robots, 1, n_timesteps
	)
	acceleration_ys = reshape(
		cat(diff(velocity_ys; dims=T), zeros(n_robots); dims=T), n_robots, 1, n_timesteps
	)
	acceleration_magnitude = reshape(
		sqrt.(acceleration_xs .^ 2 .+ acceleration_ys .^ 2), n_robots, 1, n_timesteps
	)
	angular_velocity = reshape(
		cat(
			zeros(n_robots),
			[
				abs(d) < 2pi - abs(d) ? d : -sign(d) * (2pi - abs(d)) for
				d in diff(tracking_data[:, θ, :]; dims=PROPERTIES)
			];
			dims=PROPERTIES,
		),
		n_robots,
		1,
		n_timesteps,
	)

	2 * pi .- angular_velocity[angular_velocity .> pi]
	angular_acceleration = reshape(
		cat(diff(angular_velocity; dims=T), zeros(n_robots); dims=T),
		n_robots,
		1,
		n_timesteps,
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
		dims=2,
	)

	# precalculate all metrics for all timesteps (currently no clustering)
	polarisation = [swarm_polarisation(tracking_data[:,:,t]) for t in 1:n_timesteps]
	rotational_order = [swarm_rotational_order(tracking_data[:,:,t]) for t in 1:n_timesteps]

    distmats = mapslices(
        x -> pairwise(Euclidean(), x),
        tracking_data[:, [X, Y], :]; dims=(ROBOTS, PROPERTIES),
    )
    println(typeof(distmats))
	mean_interindividual_distance = dropdims(
		mapslices(
			swarm_mean_interindividual_distance, tracking_data; dims=(ROBOTS, PROPERTIES)
		);
		dims=(ROBOTS, PROPERTIES),
	)

	surrounding_polygon =
		convex_hull.([
			[[tracking_data[r, X, t], tracking_data[r, Y, t]] for r in 1:n_robots] for
			t in 1:n_timesteps
		])
	center_of_mass = [
		[sum(tracking_data[:, X, t]), sum(tracking_data[:, Y, t])] / n_robots for
		t in 1:n_timesteps
	]
	furthest = [
		argmax(pairwise(Euclidean(), stack(surrounding_polygon[t], dims=1); dims=1)) for
		t in 1:n_timesteps
	]
	diameter = [
		Euclidean()(surrounding_polygon[t][[Tuple(furthest[t])...]]...) for
		t in 1:n_timesteps
	]
	area = [ # shoelace formula for polygon area (https://en.wikipedia.org/wiki/Shoelace_formula)
		0.5 * abs(
			sum(
				surrounding_polygon[t][i][1] *
				surrounding_polygon[t][i % length(surrounding_polygon[t]) + 1][2] -
				surrounding_polygon[t][i % length(surrounding_polygon[t]) + 1][1] *
				surrounding_polygon[t][i][2] for i in 1:length(surrounding_polygon[t])
			),
		) for t in 1:n_timesteps
	]
	roundness = [
		4 * pi * area[t] /
		sum( # perimeter
			colwise(
				Euclidean(),
				stack(surrounding_polygon[t]),
				stack(vcat(surrounding_polygon[t][2:end], [surrounding_polygon[t][1]])),
			),
		)^2 for t in 1:n_timesteps
	]

	metrics = Dict(
		"Polarisation" => polarisation,
		"Rotational Order" => rotational_order,
		"Mean IID" => mean_interindividual_distance,
		"Diameter" => diameter,
		"Area" => area,
		"Roundness" => roundness,
	)
	geometry = Dict(
		"Surrounding Polygon" => surrounding_polygon,
		"Center of Mass" => center_of_mass,
		"Furthest Robots" => furthest,
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
