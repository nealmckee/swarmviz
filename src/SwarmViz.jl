module SwarmViz

import AlgebraOfGraphics: set_aog_theme!
using Colors
using DataFrames
using Distances
using GLMakie
using JLD2
using JSON3
using LinearAlgebra
using NPZ
using NativeFileDialog
using Parquet2
using RelocatableFolders
using Statistics

# TODO julia_main function as app entry point

# constants for indexing as specified for the input data format in the readme
const AGENTS = 1
const PROPERTIES = 2
const TIME = 3
const T = 1
const X = 2
const Z = 4
const θ = 5
const HVX = 6
const HVZ = 7
const TRACKING_DIM = 5
# operating system dependent path to the assets folder
FONTFOLDER = RelocatableFolders.@path joinpath(@__DIR__, "..", "assets") #TODO make const
# color vision deficiency friendly discrete color palette
PALETTE = [ #TODO make const
	colorant"#6f4fc6",
	colorant"#ffbb00",
	colorant"#00a789",
	colorant"#d35564",
	colorant"#9799ff",
	colorant"#4b6500",
]

# load helper funtions
include("metrics.jl")
include("analysis.jl")
include("dataprocessing.jl")

"""
Struct that holds all the data input as well as everything we calculate.
"""
struct SwarmData
	"agents x properties x timesteps"
	agents::Array{Float64,3}
	"metric name => vector with metric values for each timestep"
	metrics::Dict{String,Vector{Float64}}
	"name => vector with derived data for each timestep"
	derived::Dict{String,Any}
	"one clustering object per timesteps"
	clustering::Vector{Hclust{Float64}}
end

function julia_main()::Cint
	# Set up dummy observables
	data = Observable(
		SwarmData(
			zeros(1, TRACKING_DIM + 2, 1),
			Dict("Select Metric..." => []),
			Dict(
				"Convex Hull" => [[[0.0, 0.0]]],
				"Center of Mass" => [[0.0, 0.0]],
				"Furthest Robots" => [[1, 1]],
			),
			[],
		),
	)
	wall_data = Observable(zeros(2, 1))
	wall_collisions = Observable(falses(1, 1))
	agent_collisions = Observable(falses(1, 1))
	n_timesteps = @lift size($data.agents, TIME)
	timesteps = @lift 1:($n_timesteps)
	isplaying = Observable(false)
	discrete_palette = Observable(PALETTE)

    # set the theme for the whole app
	set_aog_theme!()
	update_theme!(;
		colormap=:batlow,
		linecolor="#8f8f8f",
		markercolor="#8f8f8f",
		patchcolor="#8f8f8f",
		textcolor="#8f8f8f",
		BoxPlot=(mediancolor=:white,),
		Violin=(mediancolor=:white,),
		Figure=(backgroundcolor=:transparent,),
		Axis=(
			backgroundcolor=:transparent,
			titlecolor="#8f8f8f",
			xgridvisible=true,
			ygridvisible=true,
			xgridcolor=RGBA(143 / 255, 143 / 255, 143 / 255, 0.125),
			ygridcolor=RGBA(143 / 255, 143 / 255, 143 / 255, 0.125),
			leftspinevisible=false,
			bottomspinevisible=false,
			topspinevisible=false,
			rightspinevisible=false,
			bottomspinecolor="#8f8f8f",
			leftspinecolor="#8f8f8f",
			xtickcolor="#8f8f8f",
			ytickcolor="#8f8f8f",
			xticklabelcolor="#8f8f8f",
			yticklabelcolor="#8f8f8f",
			xlabelcolor="#8f8f8f",
			ylabelcolor="#8f8f8f",
		),
		Axis3=(
			backgroundcolor=:transparent,
			leftspinevisible=false,
			bottomspinevisible=false,
			topspinevisible=false,
			rightspinevisible=false,
			xspinecolor_1="#8f8f8f",
			yspinecolor_1="#8f8f8f",
			zspinecolor_1="#8f8f8f",
			xspinecolor_2=:transparent,
			yspinecolor_2=:transparent,
			zspinecolor_2=:transparent,
			xspinecolor_3=:transparent,
			yspinecolor_3=:transparent,
			zspinecolor_3=:transparent,
			xtickcolor="#8f8f8f",
			ytickcolor="#8f8f8f",
			ztickcolor="#8f8f8f",
			xticklabelcolor="#8f8f8f",
			yticklabelcolor="#8f8f8f",
			zticklabelcolor="#8f8f8f",
			xlabelcolor="#8f8f8f",
			ylabelcolor="#8f8f8f",
			zlabelcolor="#8f8f8f",
			titlecolor="#8f8f8f",
			xgridvisible=false,
			ygridvisible=false,
			zgridvisible=false,
			xgridcolor=RGBA(143 / 255, 143 / 255, 143 / 255, 0.125),
			ygridcolor=RGBA(143 / 255, 143 / 255, 143 / 255, 0.125),
			zgridcolor=RGBA(143 / 255, 143 / 255, 143 / 255, 0.125),
		),
		Legend=(
			framevisible=false,
			gridshalign=:left,
			padding=(0.0f0, 0.0f0, 0.0f0, 0.0f0),
			labelcolor="#8f8f8f",
			titlecolor="#8f8f8f",
		),
		Colorbar=(tickcolor="#8f8f8f",),
		size=(960, 600),
	)

    ### LAYOUT ###

	# set the title of the window that opens upon execution of the main script
	GLMakie.activate!(; title="SwarmViz")

	# create the figure and set up custom fonts for the UI
	figure = Figure(;
		size=(960, 600),
		fonts=(;
			regular=joinpath(FONTFOLDER, "fira_sans_font", "FiraSans-Medium.ttf"),
			ui_font=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Regular.ttf"),
			ui_font_bold=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Bold.ttf"),
		),
	)

	# create the four high level areas of the layout
	swarm_animation = Axis(
		figure[1, 1]; xlabel="X", ylabel="Z", autolimitaspect=1, alignmode=Outside()
	)
	metrics_grid = GridLayout(figure[1, 2]; alignmode=Outside())
	time_grid = GridLayout(figure[2, 1:2])
	controls = GridLayout(figure[3, 1:2])

	#adjust the width of the columns
	colsize!(figure.layout, 1, Relative(0.5))
	colsize!(figure.layout, 2, Relative(0.5))

	# TIMETEP

	# create the sliders that controls the current timestep together
	# with custom labels that allow changing the font
	timestep = Slider(time_grid[1, 2]; range=timesteps, startvalue=1)
	Label(time_grid[1, 1]; text="Time Step", halign=:left, font=:ui_font)
	Label(time_grid[1, 3], @lift string($(timestep.value)); halign=:right, font=:ui_font)

	# CONTROLS

	# play/pause button with a reactive label
	play_button_text = @lift $isplaying[] ? "PAUSE" : "PLAY"
	play_button = Button(
		controls[1, 1]; label=play_button_text, width=60, height=60, font=:ui_font_bold
	)

	# sliders for fps and numnber of skipped frames
	# with custom labels to allow changing the font
	animation_setting_grid = GridLayout(controls[1, 2]; default_rowgap=6)
	fps = Slider(animation_setting_grid[1, 2]; range=1:1:120, startvalue=30)
	Label(animation_setting_grid[1, 1], "FPS"; halign=:left, font=:ui_font)
	Label(
		animation_setting_grid[1, 3],
		@lift string($(fps.value));
		halign=:right,
		font=:ui_font,
	)
	skip = Slider(animation_setting_grid[2, 2]; range=0:1:99, startvalue=0)
	Label(animation_setting_grid[2, 1], "Skip"; halign=:left, font=:ui_font)
	Label(
		animation_setting_grid[2, 3],
		@lift string($(skip.value));
		halign=:right,
		font=:ui_font,
	)
	# make this column resize with the window but take up exactly 25% of the width
	colsize!(controls, 2, Relative(0.25))

	# toggles to display the convex hull, center and diameter of the swarm
	animation_toggles = [Toggle(figure) for _ in 1:3]
	controls[1, 3] = grid!(
		hcat(
			[
				Label(figure, l; halign=:left, font=:ui_font) for
				l in ["Convex Hull", "Center", "Diameter"]
			],
			animation_toggles,
		);
		default_rowgap=3,
		default_colgap=6,
		tellheight=true,
	)

	# toggles to control what is displayed on the agent markers in the animation
	# as well as which threshold is applied to the hierarchical clustering for the coloring
	agent_controls = GridLayout(controls[1, 4]; default_rowgap=3)
	agent_controls_upper = GridLayout(agent_controls[1, 1]; default_rowgap=3)
	agent_toggles = [Toggle(figure) for _ in 1:2]
	agent_controls_upper[1, 1] = grid!(
		hcat(
			[
				Label(figure, l; halign=:left, font=:ui_font) for
				l in ["Collisions", "Clustering"]
			],
			agent_toggles,
		);
		default_rowgap=3,
		default_colgap=3,
		halign=:left,
	)
	manual_log_threshold = Textbox(
		agent_controls_upper[1, 2];
		placeholder="Enter Log Threshold...",
		# font=:ui_font,
		tellwidth=false,
		halign=:right,
		reset_on_defocus=true,
		validator=Float64,
	)

	agent_controls_lower = GridLayout(agent_controls[2, 1]; default_colgap=6)
	# the possible threshold values are adapted to the minimum and maximum in the current data
	heightrange = (@lift if isempty($data.clustering)
		collect(-4:0.001:0)
	else
		round.(
			range(
				log(minimum(reduce(vcat, [c.heights for c in $data.clustering]))), 0, 1000
			),
			sigdigits=3,
		)
	end)
	log_threshold = Slider(
		agent_controls_lower[1, 2];
		range=heightrange,
		startvalue=(@lift minimum($heightrange)),
	)
	Label(agent_controls_lower[1, 1], "Log Threshold"; halign=:left, font=:ui_font)
	Label(
		agent_controls_lower[1, 3],
		@lift string($(log_threshold.value));
		halign=:right,
		font=:ui_font,
		width=25,
	)
	# make this column resize with the window and take up all remaining available space
	colsize!(controls, 4, Auto())

	# buttons for IO
	data_controls = GridLayout(controls[1, 5]; default_rowgap=6, default_colgap=6)
	import_button =
		data_controls[1:2, 1] = Button(
			figure;
			label="Import\nMovement",
			halign=:left,
			width=90,
			height=60,
			font=:ui_font,
		)
	wall_button, export_metrics_button =
		data_controls[1:2, 2] = [
			Button(figure; label=l, halign=:left, width=w, height=27, font=:ui_font) for
			(l, w) in zip(["Import Wall", "Export All"], [90, 90])
		]

	# METRIC PLOTS

	# axes to hold the metric plots with a narrow one below them for the collision events
	metric_axes = [Axis(metrics_grid[row, 1]; xticklabelsvisible=false) for row in 1:3]
	collisions_axis = Axis(
		metrics_grid[4, 1];
		xlabel="Time Step",
		yticklabelsvisible=false,
		yticksvisible=false,
		height=10,
		ygridvisible=false,
		xgridvisible=false,
	)
	linkxaxes!(metric_axes..., collisions_axis)

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

    ### INTERACTIVITY ###

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
		JSON3.write(
			joinpath(export_folder, "convex_hull.json"), data[].derived["Convex Hull"]
		)
	end

    ### PLOTTING ###

	# ANIMATION

	# Make coordinates and rotation responsive to the time slider
	x, y, r = [@lift $data.agents[:, i, $(timestep.value)] for i in [X, Z, θ]]
	# make color of agents dependent on clustering
	agent_colors = cluster_coloring_obs(
		data, timestep, log_threshold, agent_toggles, discrete_palette
	)
	# make glowcolor of agents dependend on collisions
	g = collision_coloring_obs(
		data, agent_collisions, wall_collisions, timestep, skip, agent_toggles
	)

	# Plot the agent swarm
	agent_marker = Makie.Polygon(Point2f[(-4, -3), (-1, 0), (-4, 3), (5, 0)])
	marker_scale = (@lift round(Int, 6 / sqrt(size($data.agents, 1))))
	glow_scale = (@lift $marker_scale * 3)
	scatter!(
		swarm_animation,
		x,
		y;
		marker=agent_marker,
		markersize=marker_scale,
		rotations=r,
		color=agent_colors,
		glowcolor=g,
		glowwidth=glow_scale,
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
				$data.agents[
					$data.derived["Furthest Robots"][$(timestep.value)],
					[X, Z],
					$(timestep.value),
				];
				dims=AGENTS,
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
			width=190,
			height=24,
			halign=:left,
			selection_cell_color_inactive=:transparent,
			textcolor="#8f8f8f",
			prompt="Select Metric...",
		) for i in 1:3
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
		connect!(p.visible, agent_toggles[1].active)
	end

	# plot timestep markers in metrics and collision axes
	for axis in metric_axes
		vlines!(axis, timestep.value; color=:black, linewidth=0.5)
	end

	# Display the resulting figure in its own window
	display(figure)
	return 0
end

end
