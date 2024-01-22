function swarm_polarisation(tracking_datapoint)
	return norm(sum(tracking_datapoint[:, [HVX, HVZ]]; dims=1)) /
		   size(tracking_datapoint, 1)
end

# Calculate the rotational order of the swarm for the data of one timestep
function swarm_rotational_order(tracking_datapoint)
	clustersize = size(tracking_datapoint, ROBOTS)
	barycenter = sum(tracking_datapoint[:, [X, Z]]; dims=ROBOTS) ./ clustersize
	radial_unit_vectors = mapslices(
		x -> x / norm(x), tracking_datapoint[:, [X, Z]] .- barycenter; dims=PROPERTIES
	)
	return norm(
		sum([
			cross(
				vcat(tracking_datapoint[i, [HVX, HVZ]], 0),
				vcat(radial_unit_vectors[i, :], 0),
			) for i in 1:clustersize
		]),
	) / clustersize
end

function swarm_mean_interindividual_distance(tracking_datapoint)
	return mean(pairwise(Euclidean(), tracking_datapoint[:, [X, Z]]; dims=ROBOTS))
end #TODO reuse Euclidean distance for futhest, diameter, mean lowest

function swarm_max_interindividual_distance(tracking_datapoint) #TODO unused -> remove
	return maximum(pairwise(Euclidean(), tracking_datapoint[:, [X, Z]]; dims=ROBOTS))
end

#TODO: add mean lowest distance?
