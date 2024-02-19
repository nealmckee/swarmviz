function swarm_polarization(tracking_datapoint)
	return norm(sum(tracking_datapoint[:, [HVX, HVZ]]; dims=1)) /
		   size(tracking_datapoint, 1)
end

# Calculate the rotational order of the swarm for the data of one timestep
function swarm_rotational_order(tracking_datapoint, center_of_mass)
	clustersize = size(tracking_datapoint, AGENTS)
	radial_unit_vectors = mapslices(
		x -> x / norm(x), tracking_datapoint[:, [X, Z]] .- center_of_mass'; dims=PROPERTIES
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

function shoelace_area(polygon)
	return 0.5 * abs(
		sum(a[1] * b[2] - b[1] * a[2] for (a, b) in zip(polygon, circshift(polygon, -1)))
	)
end
