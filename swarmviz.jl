using LinearAlgebra
using GLMakie
using NPZ
using Distances
using Statistics

experiment_file = "data/DatasetE2/E21/E212/E212r1_summaryd.npy"

# Load the data and drop singular dimension of runs
tracking_data = npzread(experiment_file)[1, :, :, :]

# Get the number of robots and the number of timesteps
n_robots = size(tracking_data, 1)
n_timesteps = size(tracking_data, 3)

# Add heading vector from orientation to data
tracking_data = cat(
    tracking_data, reshape(cos.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps); dims=2
)
tracking_data = cat(
    tracking_data, reshape(sin.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps); dims=2
)

# Calculate the polarisation of the swarm for the data of one timestep
function swarm_polarisation(tracking_datapoint; cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    return norm(sum(tracking_datapoint[:, 6:7]; dims=1)) / size(tracking_datapoint, 1)
end

# Calculate the rotational order of the swarm for the data of one timestep
function swarm_rotational_order(tracking_datapoint; cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    clustersize = size(tracking_datapoint, 1)
    barycenter = sum(tracking_datapoint[:, 1:2]; dims=1) ./ clustersize
    radial_unit_vectors = mapslices(
        x -> x / norm(x), tracking_datapoint[:, 1:2] .- barycenter; dims=2
    )
    return norm(
        sum([
            cross(
                push!(tracking_datapoint[i, 6:7], 0), push!(radial_unit_vectors[i, :], 0)
            ) for i in 1:clustersize
        ]),
    ) / clustersize
end

function swarm_mean_interindividual_distance(tracking_datapoint, cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    return mean(pairwise(Euclidean(), tracking_datapoint[:, 1:2]; dims=2))
end

# precalculate all metrics for all timesteps (currently no clustering)
polarisation = dropdims(
    mapslices(swarm_polarisation, tracking_data; dims=(1, 2)); dims=(1, 2)
)
rotational_order = dropdims(
    mapslices(swarm_rotational_order, tracking_data; dims=(1, 2)); dims=(1, 2)
)
mean_interindividual_distance = dropdims(
    mapslices(swarm_mean_interindividual_distance, tracking_data; dims=(1, 2)); dims=(1, 2)
)
