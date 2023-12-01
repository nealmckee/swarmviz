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
    barycenter = sum(tracking_datapoint[:, [2, 4]]; dims=1) ./ clustersize
    radial_unit_vectors = mapslices(
        x -> x / norm(x), tracking_datapoint[:, [2, 4]] .- barycenter; dims=2
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
    return mean(pairwise(Euclidean(), tracking_datapoint[:, [2, 4]]; dims=1))
end

function swarm_max_interindividual_distance(tracking_datapoint, cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    return maximum(pairwise(Euclidean(), tracking_datapoint[:, [2, 4]]; dims=1))
end
