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

# set the log threshold to the content of the textbox if it can be parsed to a Float64
on(manual_log_threshold.stored_string) do s
	set_close_to!(log_threshold, parse(Float64, s))
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
			size(data[].agents, 1) - length(PALETTE),
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
		minimum(data[].agents[:, X, :]) - 100,
		maximum(data[].agents[:, X, :]) + 100,
		minimum(data[].agents[:, Z, :]) - 100,
		maximum(data[].agents[:, Z, :]) + 100,
	)
	# load and process the collision data if it’s present in the same directory
	for (name, obs) in
		zip(["warefl.json", "aarefl.json"], [wall_collisions, agent_collisions])
		path = replace(tracking_path, "summaryd.npy" => name)
		if isfile(path)
			obs[] = process_collisions(path)
		end
	end
	# adjust the limits of the collisions plot (and with that also the metrics plots)
	xlims!(collisions_axis, 1, n_timesteps[])
end

on(wall_button.clicks) do c
	# on button press, opens a native file dialogue
	wall_path = pick_file(; filterlist="npy")
	# then, if the dialog wasn’t cancelled, loads and analyses the chosen .npy file
	wall_path == "" || (wall_data[] = analyse_wall(wall_path))
	# adjust the limits of the animation to fit the new wall data
	autolimits!(swarm_animation)
end

on(export_metrics_button.clicks) do c
	# on button press, opens a native file dialogue
	export_folder = pick_folder()
	# if it’s not cancelled
	export_folder == "" && return nothing
	# first reshape the datapoints into their different output formats if necessary
	agents_df = agentdata2longerdf(data, log_threshold)
	tensors = derived2tensors(data)
	chosen_clustering_threshold = exp(log_threshold.value[])
	# then write them to the matching files types
	Parquet2.writefile(joinpath(export_folder, "agents.parquet"), agents_df)
	Parquet2.writefile(
		joinpath(export_folder, "metrics.parquet"), DataFrame(data[].metrics)
	)
	jldopen(joinpath(export_folder, "derived.jld2"), "w") do file
		for (p, t) in tensors
			file[p] = t
		end
		file["clustering/chosen_clustering_threshold"] = chosen_clustering_threshold
	end
	JSON3.write(joinpath(export_folder, "convex_hull.json"), data[].derived["Convex Hull"])
end
