# This scripts set ups all the reactive elements, mainly the plots themselves

# ANIMATION

# Make coordinates and rotation responsive to the time slider
x, y, r = [@lift $data.robots[:, i, $(timestep.value)] for i in [X, Z, θ]]
# make color of robots dependent on clustering
robot_colors = cluster_coloring_obs(
	data, timestep, log_threshold, robot_toggles, discrete_palette
)
# make glowcolor of robots dependend on collisions
g = collision_coloring_obs(
	data, agent_collisions, wall_collisions, timestep, skip, robot_toggles
)

# Plot the robot swarm
robot_marker = Makie.Polygon(Point2f[(-4, -3), (-1, 0), (-4, 3), (5, 0)])
marker_scale = (@lift round(Int, 6 / sqrt(size($data.robots, 1))))
scatter!(
	swarm_animation,
	x,
	y;
	marker=robot_marker,
	markersize=marker_scale,
	rotations=r,
	color=robot_colors,
	glowcolor=g,
	glowwidth=marker_scale,
)

# Plot the wall of the enclosure
wall_vertices = @lift Point2f.($wall_data[1, :], $wall_data[2, :])
poly!(
	swarm_animation,
	wall_vertices;
	color=:transparent,
	strokecolor="#8f8f8f",
	strokewidth=1,
	linestyle=:dot,
	closed=true,
	depth=-1,
)

# plot the surrounding polygon and connect to toggle
surrounding_polygon = poly!(
	swarm_animation,
	@lift Point2f.($data.derived["Convex Hull"][$(timestep.value)]);
	color=:transparent,
	strokecolor="#8f8f8f",
	strokewidth=1,
	linestyle=:dash,
	closed=true,
)
connect!(surrounding_polygon.visible, animation_toggles[1].active)

# plot the center of mass and connect to toggle
center_of_mass = scatter!(
	swarm_animation,
	@lift Point2f($data.derived["Center of Mass"][$(timestep.value)]);
	color="#8f8f8f",
	markersize=12,
	marker=:xcross,
)
connect!(center_of_mass.visible, animation_toggles[2].active)

#plot the diameter and connect to toggle
diameter = lines!(
	swarm_animation,
	@lift Point2f.(
		eachslice(
			$data.robots[
				$data.derived["Furthest Robots"][$(timestep.value)],
				[X, Z],
				$(timestep.value),
			];
			dims=ROBOTS,
		)
	);
	color="#8f8f8f",
	linewidth=1,
)
connect!(diameter.visible, animation_toggles[3].active)

# METRICS

# menus to select which metrics to plot (inside the metric axes)
metric_tuples = @lift Tuple.(collect($data.metrics))
metric_menus = [
	Menu(
		metrics_grid[i, 1, Top()];
		options=metric_tuples,
		default=default,
		width=150, #TODO: adapt length to max length of final metric selection
		height=24,
		halign=:left,
		selection_cell_color_inactive=:transparent,
		textcolor="#8f8f8f",
		prompt="Select Metric...",
	) for (i, default) in enumerate(["Polarisation", "Rotational Order", "Diameter"]) #TODO remove defaults after debugging
]

# plot the metrics when one is chosen from the corresponding menu
for (i, (menu, axis)) in enumerate(zip(metric_menus, metric_axes))
	on(menu.selection) do s
		while length(axis.scene.plots) > 1
			delete!(
				axis,
				axis.scene.plots[typeof.(axis.scene.plots) .!= Plot{
					Makie.vlines,Tuple{Int64}
				}][1],
			)
		end
		if !isnothing(s)
			lines!(axis, s; linewidth=1, color=PALETTE[i])
			limits!(
				axis,
				(nothing, nothing),
				maximum(s) <= 1 && minimum(s) >= 0 ? (0, 1) : (nothing, nothing),
			)
		end
	end
	notify(menu.selection) #TODO remove after debugging
end

# plot collisions
collision_plots = [
	vlines!(
		collisions_axis,
		@lift findall(reduce(|, $obs; dims=1)[1, :]);
		color=color,
		linewidth=1,
		depth=depth,
	) for (depth, (obs, color)) in
	enumerate(zip([wall_collisions, agent_collisions], [PALETTE[4], PALETTE[5]]))
]
for p in collision_plots
	connect!(p.visible, robot_toggles[1].active)
end

# plot timestep markers in metrics and collision axes
for axis in metric_axes
	vlines!(axis, timestep.value; color=:black, linewidth=0.5)
end
