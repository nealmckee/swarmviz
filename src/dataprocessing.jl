function agentdata2longerdf(data, log_threshold)
	df = DataFrame(
		reshape(permutedims(data[].agents, (1, 3, 2)), (:, size(data[].agents, 2), 1))[
			:, :, 1
		],
		[:t, :x, :y, :z, :angle_xz, :heading_x, :heading_z],
	)
	df.agent_id = repeat(1:size(data[].agents, 1); inner=size(data[].agents, 3))
	df[!, "Clustering at Threshold = $(exp(log_threshold.value[]))"] = reduce(
		vcat, cutree.(data[].clustering; h=log_threshold.value[])
	)
	return df
end

function derived2tensors(data)
	tensors = Dict()
	tensors["distance_matrices"] = stack(data[].derived["Distance Matrices"])
	tensors["furthest_agents"] = stack(data[].derived["Furthest Robots"])
	tensors["center_of_mass"] = stack(data[].derived["Center of Mass"])
	tensors["clustering/clustering_heights"] = stack([c.heights for c in data[].clustering])
	tensors["clustering/clustering_merges"] = stack([c.merges for c in data[].clustering])
	tensors["clustering/clustering_order"] = stack([c.order for c in data[].clustering])
	return tensors
end

function cluster_coloring_obs(
	data, timestep, log_threshold, animation_toggles, discrete_palette
)
	return @lift (
		if !$(animation_toggles[6].active) || data[].clustering == []
			repeat([RGBA(0, 0, 0, 1)], size($data.agents, 1))
		else
			$discrete_palette[collect(
				cutree($data.clustering[$(timestep.value)]; h=exp($(log_threshold.value)))
			)]
		end
	)
end

function recent_collision(collisions, agent_index, timestep_value, skip_value)
	!checkbounds(Bool, collisions, 1, timestep_value) && return false
	collision_occured = any(
		collisions[agent_index, max((timestep_value - skip_value), 1):(timestep_value)]
	)
	return collision_occured
end

function collision_coloring_obs(
	data, agent_collisions, wall_collisions, timestep, skip, animation_toggles
)
	return @lift (
		if size($agent_collisions, 1) != size($data.agents, 1) ||
			size($wall_collisions, 1) != size($data.agents, 1) ||
			!$(animation_toggles[4].active)
			repeat([RGBA(0, 0, 0, 0)], size($data.agents, 1))
		else
			[
				if recent_collision($agent_collisions, i, $(timestep.value), $(skip.value))
					PALETTE[5]
				elseif recent_collision(
					$wall_collisions, i, $(timestep.value), $(skip.value)
				)
					PALETTE[4]
				else
					RGBA(0, 0, 0, 0)
				end for i in axes($wall_collisions, 1)
			]
		end
	)
end
