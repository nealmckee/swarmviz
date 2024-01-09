function analyse_tracking(filename)
	# Load the data and drop singular dimensions
	tracking_data = npzread(filename)[1, :, 1:5, :]

	# Get the number of robots and the number of timesteps
	n_robots = size(tracking_data, 1)
	n_timesteps = size(tracking_data, 3)

	# Add heading vector from orientation to data
	heading_vector_xs = reshape(cos.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps)
	heading_vector_ys = reshape(sin.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps)
	velocity_xs = reshape(
		cat(zeros(10), diff(tracking_data[:, 2, :]; dims=2); dims=2),
		n_robots,
		1,
		n_timesteps,
	)
	velocity_ys = reshape(
		cat(zeros(10), diff(tracking_data[:, 2, :]; dims=2); dims=2),
		n_robots,
		1,
		n_timesteps,
	)
	velocity_magnitude = reshape(
		sqrt.(velocity_xs .^ 2 .+ velocity_ys .^ 2), n_robots, 1, n_timesteps
	)
	acceleration_xs = reshape(
		cat(diff(velocity_xs; dims=2), zeros(10); dims=2), n_robots, 1, n_timesteps
	)
	acceleration_ys = reshape(
		cat(diff(tracking_data[:, 2, :]; dims=2), zeros(10); dims=2),
		n_robots,
		1,
		n_timesteps,
	)
	acceleration_magnitude = reshape(
		sqrt.(acceleration_xs .^ 2 .+ acceleration_ys .^ 2), n_robots, 1, n_timesteps
	)
	angular_velocity = reshape(
		cat(zeros(10), diff(tracking_data[:, 5, :]; dims=2); dims=2),
		n_robots,
		1,
		n_timesteps,
	)
	angular_acceleration = reshape(
		cat(diff(angular_velocity; dims=2), zeros(10); dims=2), n_robots, 1, n_timesteps
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
	polarisation = dropdims(
		mapslices(swarm_polarisation, tracking_data; dims=(1, 2)); dims=(1, 2)
	)
	rotational_order = dropdims(
		mapslices(swarm_rotational_order, tracking_data; dims=(1, 2)); dims=(1, 2)
	)
	mean_interindividual_distance = dropdims(
		mapslices(swarm_mean_interindividual_distance, tracking_data; dims=(1, 2));
		dims=(1, 2),
	)
	metrics = transpose(
		cat(polarisation, rotational_order, mean_interindividual_distance; dims=2)
	)
	return SwarmData(tracking_data[:, 1:5, :], tracking_data[:, 5:end, :], metrics)
end

function analyse_wall(filename)
	wd = npzread(filename)[1, 1, [2, 4], :]
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
