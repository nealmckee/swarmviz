function analyse_tracking(filename)
    # Load the data and drop singular dimensions
    tracking_data = npzread(filename)[1, :, :, :]

    # Get the number of robots and the number of timesteps
    n_robots = size(tracking_data, 1)
    n_timesteps = size(tracking_data, 3)

    # Add heading vector from orientation to data
    tracking_data = cat(
        tracking_data,
        reshape(cos.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps);
        dims=2,
    )
    tracking_data = cat(
        tracking_data,
        reshape(sin.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps);
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
    metrics = cat(
        polarisation,
        rotational_order,
        mean_interindividual_distance;
        dims=2,
    ) |> transpose
    return SwarmData(tracking_data, metrics)
end

function analyse_wall(filename)
    wd = npzread(filename)[1, 1, [2, 4], :]
    # Remove columns with NaN values from wall_data
    wd = wd[:, .!isnan.(wd[1, :])]
    wd = wd[:, .!isnan.(wd[2, :])]
    return wd
end
