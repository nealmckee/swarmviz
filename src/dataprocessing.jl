function robotdata2longerdf(data, threshold)
	df = DataFrame(
		reshape(permutedims(data[].robots, (1, 3, 2)), (:, size(data[].robots, 2), 1))[
			:, :, 1
		],
		[ #TODO check what’s left after discussing preprocessing etc.
			:t,
			:x,
			:y,
			:z,
			:angle_xz,
			:heading_x,
			:heading_z,
			:velocity_x,
			:velocity_z,
			:l2norm_velocity_xy,
			:acceleration_x,
			:acceleration_z,
			:l2norm_acceleration_xz,
			:angular_velocity,
			:angular_acceleration,
		],
	)
	df.robot_id = repeat(1:size(data[].robots, 1); inner=size(data[].robots, 3))
	df[!, "Clustering at Threshold = $(threshdold.value[])"] = reduce(
		vcat, cutree.(data[].clustering; h=threshold.value[])
	)
	return df
end

function derived2tensors(data)
	distance_matrices = stack(data[].derived["Distance Matrices"])
	furthest_robots = stack(data[].derived["Furthest Robots"])
	center_of_mass = stack(data[].derived["Center of Mass"])
	return distance_matrices, furthest_robots, center_of_mass
end

function cluster_coloring_obs(data, timestep, threshold, robot_toggles, discrete_palette)
	return @lift (
		if !$(robot_toggles[2].active)
			repeat([RGBA(0, 0, 0, 1)], size($data.robots, 1))
		else
			$discrete_palette[collect(
				cutree($data.clustering[$(timestep.value)]; h=$(threshold.value))
			)]
		end
	)
end

function collision_coloring_obs(
	data, agent_collisions, wall_collisions, timestep, skip, robot_toggles
)
	return @lift ( #TODO: refactor?
		if size($agent_collisions, 1) != size($data.robots, 1) ||
			size($wall_collisions, 1) != size($data.robots, 1) ||
			!$(robot_toggles[1].active)
			repeat([RGBA(0, 0, 0, 0)], size($data.robots, 1))
		else
			[
				if checkbounds(Bool, $agent_collisions, 1, $(timestep.value)) && any(
					$agent_collisions[
						i, max(($(timestep.value) - skip.value[]), 1):($(timestep.value))
					],
				)
					PALETTE[5]
				elseif checkbounds(Bool, $wall_collisions, 1, $(timestep.value)) && any(
					$wall_collisions[
						i, max(($(timestep.value) - skip.value[]), 1):($(timestep.value))
					],
				)
					PALETTE[4]
				else
					RGBA(0, 0, 0, 0)
				end for i in 1:size($wall_collisions, 1)
			]
		end
	)
end
