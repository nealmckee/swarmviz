function analyse_tracking(filename)
    # Load the data and drop singular dimensions
    td = npzread(filename)[1, :, :, :]

    # Get the number of robots and the number of timesteps
    n_robots = size(td, 1)
    nt = size(td, 3)

    # Add heading vector from orientation to data
    td = cat(td, reshape(cos.(td[:, 5, :]), n_robots, 1, nt); dims=2)
    td = cat(td, reshape(sin.(td[:, 5, :]), n_robots, 1, nt); dims=2)

    # precalculate all metrics for all timesteps (currently no clustering)
    pol = dropdims(mapslices(swarm_polarisation, td; dims=(1, 2)); dims=(1, 2))
    ro = dropdims(
        mapslices(swarm_rotational_order, td; dims=(1, 2)); dims=(1, 2)
    )
    miid = dropdims(
        mapslices(swarm_mean_interindividual_distance, td; dims=(1, 2)); dims=(1, 2)
    )
    return td, pol, ro, miid, nt
end

function analyse_wall(filename)
    wd = npzread(filename)[1, 1, [2, 4], :]
    # Remove columns with NaN values from wall_data
    wd = wd[:, .!isnan.(wd[1, :])]
    wd = wd[:, .!isnan.(wd[2, :])]
    return wd
end
