# This script adds interactivity to all the previously created elements of the layout.

# TIMESTEP

# starts animation loop on buttonpress, the plot then automatically updates
# as it’s dependent on the time slider
on(play_button.clicks) do c
	isplaying[] = !isplaying[]
end
on(play_button.clicks) do c
	@async while isplaying[] && timestep.value[] < n_timesteps[] - skip.value[] - 1
		frame_start = time()
		set_close_to!(timestep, timestep.value[] + 1 + skip.value[])
		sleep(max(1 / fps.value[] - time() + frame_start, 0))
		isopen(figure.scene) || break # ensures computations stop if window is closed
		if timestep.value[] >= n_timesteps[] - skip.value[] - 1
			isplaying[] = !isplaying[]
		end
	end
end

# CONTROLS

# set the threshold to the content of the textbox if it can be parsed to a Float64
on(manual_threshold.stored_string) do s
	set_close_to!(threshold, parse(Float64, s))
end

# IO

on(import_button.clicks) do c
	# on button press, opens a native file dialogue
	tracking_path = pick_file(; filterlist="npy")
	# preemptively set the timestep to 1
	# as the plots will be redrawn immediately after changing the data observable
	set_close_to!(timestep, 1)
	# then, if the dialog wasn’t cancelled, loads and analyses the chosen .npy file
	tracking_path == "" && return nothing
	data[] = analyse_tracking(tracking_path)
	discrete_palette[] = vcat(
		RGB.(PALETTE),
		distinguishable_colors(
			size(data[].robots, 1) - length(PALETTE),
			vcat(RGB.(PALETTE), [RGB(0, 0, 0), RGB(1, 1, 1)]);
			dropseed=true,
			lchoices=range(45, 90; length=15),
		),
	)
	# adjust the limits of the plots to the new data
	autolimits!(metric_axes[1])
	# (with padding for the swarm animation)
	limits!(
		swarm_animation,
		minimum(data[].robots[:, X, :]) - 100,
		maximum(data[].robots[:, X, :]) + 100,
		minimum(data[].robots[:, Z, :]) - 100,
		maximum(data[].robots[:, Z, :]) + 100,
	)
end

on(wall_button.clicks) do c
	# on button press, opens a native file dialogue
	wall_path = pick_file(; filterlist="npy")
	# then, if the dialog wasn’t cancelled, loads and analyses the chosen .npy file
	wall_path == "" || (wall_data[] = analyse_wall(wall_path))
	# adjust the limits of the animation to fit the new wall data
	autolimits!(swarm_animation)
end

on(collision_button.clicks) do c
	# on button press, opens a native file dialogue
	collisions_folder = pick_folder()
	# parses the chosen directory for files with the ending specified in the readme
	wall_collisons_path = filter(
		x -> occursin(r".*warefl.json", x), readdir(collisions_folder; join=true)
	)[1]
	agent_collisons_path = filter(
		x -> occursin(r".*aarefl.json", x), readdir(collisions_folder; join=true)
	)[1]
	# load and process the collision data
	wall_collisions[] = process_collisions(wall_collisons_path)
	agent_collisions[] = process_collisions(agent_collisons_path)
	# adjust the limits of the plot to fit the new data
	# (if it’s currently displayed)
	autolimits!(collisions_axis)
end

on(export_metrics_button.clicks) do c
	# on button press, opens a native file dialogue
	export_folder = pick_folder()
	# if it’s not cancelled
	export_folder != "" && return nothing
	# first reshape the datapoints into their different output formats if necessary
	robots_df = robotdata2longerdf(data, threshold)
	distance_matrices, furthest_robots, center_of_mass = derived2tensors(data)
	chosen_clustering_threshold = threshold.value[]
	# then write them to the matching files types
	Parquet2.writefile(joinpath(export_folder, "robots.parquet"), robots_df)
	Parquet2.writefile(
		joinpath(export_folder, "metrics.parquet"), DataFrame(data[].metrics)
	)
	JLD2.@save joinpath(export_folder, "derived.jld2") distance_matrices furthest_robots center_of_mass chosen_clustering_threshold
	JSON3.write(joinpath(export_folder, "convex_hull.json"), data[].derived["Convex Hull"])
end
