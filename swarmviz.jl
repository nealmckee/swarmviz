import AlgebraOfGraphics: set_aog_theme!
using Colors
using DataFrames
using Distances
using GLMakie
using JSON3
using JLD2
using LinearAlgebra
using NativeFileDialog
using NPZ
using Parquet2
using RelocatableFolders
using Statistics

const ROBOTS = 1
const PROPERTIES = 2
const T = 3
const X = 2
const Z = 4
const θ = 5
const HVX = 6
const HVZ = 7
const ACCX = 11
const ACCZ = 12
const TRACKING_DIM = 5
const FONTFOLDER = RelocatableFolders.@path joinpath(@__DIR__, "assets")

include("src/metrics.jl")
include("src/plotstyle.jl")
include("src/analysis.jl")

"""
Struct that holds all the data input as well as everything we calculate.
"""
struct SwarmData
	"robots x properties x timesteps"
	robots::Array{Float64,3}
	"metric name => vector with metric values for each timestep"
	metrics::Dict{String,Union{Vector{Float64},Matrix{Float64}}}
	"name => vector with derived ata for each timestep"
	derived::Dict{String,Any}
	"one clustering per timesteps"
	clustering::Vector{Hclust{Float64}}
end

# Set up observables
#TODO: check without dummy data
data = Observable(SwarmData(zeros(1, TRACKING_DIM + 10, 1), Dict(), Dict(), []))
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(falses(1, 1))
agent_collisions = Observable(falses(1, 1))
n_timesteps = @lift size($data.robots, T)
timesteps = @lift 1:($n_timesteps)
isplaying = Observable(false)
discrete_palette = Observable(Makie.wong_colors())

# Preload data for easier debugging #TODO: remove when done
tracking_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_summaryd.npy"
wall_path = "data/A0_09_B0_0/w3_summaryd.npy"
wall_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_warefl.json"
agent_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_aarefl.json"
wall_data[] = analyse_wall(wall_path)
data[] = analyse_tracking(tracking_path)
wall_collisions[] = process_collisions(wall_collisons_path)
agent_collisions[] = process_collisions(agent_collisons_path)
discrete_palette[] = distinguishable_colors(
	size(data[].robots, 1), RGB.(Makie.wong_colors())
)

# Set up the figure
GLMakie.activate!(; title="SwarmViz")
fig = Figure(;
	size=(960, 600),
	fonts=(;
		regular=joinpath(FONTFOLDER, "fira_sans_font", "FiraSans-Medium.ttf"),
		ui_font=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Regular.ttf"),
		ui_font_bold=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Bold.ttf"),
	),
)

swarm_animation = Axis(
	fig[1, 1]; xlabel="X", ylabel="Z", autolimitaspect=1, alignmode=Outside()
)
metrics_grid = GridLayout(fig[1, 2]; alignmode=Outside())
time_grid = GridLayout(fig[2, 1:2])
controls = GridLayout(fig[3, 1:2])

#adjust the width of the column with the swarm animation
colsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))

Label(time_grid[1, 1]; text="Timestep", halign=:left, font=:ui_font)
timestep = Slider(time_grid[1, 2]; range=timesteps, startvalue=1)
Label(time_grid[1, 3], @lift string($(timestep.value)); halign=:right, font=:ui_font)
# TODO: test text input for time
# on(timestep_text.stored_string) do s
# 	timestep.value[] = parse(Int64, s)
# end
# on(timestep.value) do v
# 	timestep_text.stored_string.val = string(v)
# end

# create play/pause button with reactive label
play_button_text = @lift $isplaying[] ? "PAUSE" : "PLAY"
play_button = Button(
	controls[1, 1]; label=play_button_text, width=60, height=60, font=:ui_font_bold
)

animation_setting_grid = GridLayout(controls[1, 2]; default_rowgap=6)
fps = Slider(animation_setting_grid[1, 2]; range=1:1:120, startvalue=30)
Label(animation_setting_grid[1, 1], "FPS"; halign=:left, font=:ui_font)
Label(
	animation_setting_grid[1, 3], @lift string($(fps.value)); halign=:right, font=:ui_font
)
skip = Slider(animation_setting_grid[2, 2]; range=0:1:99, startvalue=0)
Label(animation_setting_grid[2, 1], "Skip"; halign=:left, font=:ui_font)
Label(
	animation_setting_grid[2, 3], @lift string($(skip.value)); halign=:right, font=:ui_font
)
colsize!(controls, 2, Relative(0.25))

animation_toggles = [Toggle(fig) for _ in 1:3]
controls[1, 3] = grid!(
	hcat(
		[
			Label(fig, l; halign=:left, font=:ui_font) for
			l in ["Convex Hull", "Center", "Diameter"]
		],
		animation_toggles,
	);
	default_rowgap=3,
	default_colgap=6,
	tellheight=true,
)

robot_controls = GridLayout(controls[1, 4]; default_rowgap=3)
robot_toggles = [Toggle(fig) for _ in 1:2]
robot_controls[1, 1:3] = grid!(
	hcat(
		[Label(fig, l; halign=:left, font=:ui_font) for l in ["Collisions", "Clustering"]],
		robot_toggles,
	);
	default_rowgap=3,
	default_colgap=15,
	halign=:left,
)
heightrange = @lift round.(
	range(extrema(reduce(vcat, [c.heights for c in $data.clustering]))..., 3000),
	sigdigits=4,
)
threshold = Slider(
	robot_controls[2, 2]; range=heightrange, startvalue=(@lift median($heightrange))
)
Label(robot_controls[2, 1], "Threshold"; halign=:left, font=:ui_font)
Label(robot_controls[2, 3], @lift string($(threshold.value)); halign=:right, font=:ui_font)
colsize!(controls, 4, Auto())

# on(threshold.value) do v # TODO remove after feedback
# 	_mean_clustersize = [
# 		mean(values(countmap(cutree(c; h=v)))) for c in data[].clustering
# 	]
# 	# idcs = findall(m -> m.i_selected[] == 8, metric_menus)
# 	# for (i, menu, axis) in zip(idcs, metric_menus[idcs], metric_axes[idcs])
# 	# 	empty!(axis)
#     #     lines!(axis, _mean_clustersize; linewidth=1, color=Makie.wong_colors()[i])
#     #     autolimits!(axis)
#     #     vlines!(axis, timestep.value; color=:black, linewidth=0.5)
# 	# end
# end

# start animation loop on buttonpress, the plot then automatically updates
# as it’s dependent on the time slider
on(play_button.clicks) do c
	isplaying[] = !isplaying[]
end
on(play_button.clicks) do c
	@async while isplaying[] && timestep.value[] < n_timesteps[] - skip.value[] - 1
		frame_start = time()
		set_close_to!(timestep, timestep.value[] + 1 + skip.value[])
		sleep(max(1 / fps.value[] - time() + frame_start, 0))
		isopen(fig.scene) || break # ensures computations stop if window is closed
		if timestep.value[] >= n_timesteps[] - skip.value[] - 1
			isplaying[] = !isplaying[]
		end
	end
end

# create buttons to import/export data
data_controls = GridLayout(controls[1, 5]; default_rowgap=6, default_colgap=6)
import_button, wall_button, collision_button, export_metrics_button =
	data_controls[1:2, 1:2] = [
		Button(fig; label=l, halign=:left, width=w, height=27, font=:ui_font) for
		(l, w) in zip(
			["Import Tracking", "Import Wall", "Import Collisions", "Export All"],
			[120, 90, 120, 90],
		)
	]

on(import_button.clicks) do c
	tracking_path = pick_file(; filterlist="npy")
	if tracking_path != ""
		data[] = analyse_tracking(tracking_path)
		discrete_palette[] = distinguishable_colors( # TODO: add black and white to seed
			size(data[].robots, 1),
			RGB.(Makie.wong_colors()),
		)
	end
	autolimits!(metric_axes[1])
	limits!(
		swarm_animation,
		minimum(data[].robots[:, X, :]) - 100,
		maximum(data[].robots[:, X, :]) + 100,
		minimum(data[].robots[:, Z, :]) - 100,
		maximum(data[].robots[:, Z, :]) + 100,
	)
	set_close_to!(timestep, 1)
end

on(wall_button.clicks) do c
	wall_path = pick_file(; filterlist="npy")
	wall_path == "" || (wall_data[] = analyse_wall(wall_path))
	autolimits!(swarm_animation)
end

on(collision_button.clicks) do c
	collisions_folder = pick_folder()
	wall_collisons_path = filter(
		x -> occursin(r".*warefl.json", x), readdir(collisions_folder; join=true)
	)[1]
	agent_collisons_path = filter(
		x -> occursin(r".*aarefl.json", x), readdir(collisions_folder; join=true)
	)[1]
	wall_collisions[] = process_collisions(wall_collisons_path)
	agent_collisions[] = process_collisions(agent_collisons_path)
	autolimits!(collisions_axis)
end

on(export_metrics_button.clicks) do c
	export_folder = pick_folder()
	robots_df = DataFrame(
		reshape(permutedims(data[].robots, (1, 3, 2)), (:, size(data[].robots, 2), 1))[
			:, :, 1
		],
		[
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
	robots_df.robot_id = repeat(1:size(data[].robots, 1); inner=size(data[].robots, 3))
	Parquet2.writefile("robots.parquet", robots_df)
	Parquet2.writefile("metrics.parquet", DataFrame(data[].metrics))
	distance_matrices = stack(data[].derived["Distance Matrices"])
	furthest_robots = stack(data[].derived["Furthest Robots"])
	center_of_mass = stack(data[].derived["Center of Mass"])
	JLD2.@save "derived.jld2" distance_matrices furthest_robots center_of_mass
	JSON3.write("convex_hull.json", data[].derived["Convex Hull"])
end

# axes to hold the metric plots
metric_axes = [
	Axis(
		metrics_grid[row, 1];
		xticklabelsvisible=false,
		palette=(patchcolor=collect(cgrad(:batlow, 9; categorical=true)),),
	) for row in 1:3
]
collisions_axis = Axis(
	metrics_grid[4, 1];
	xlabel="Timestep",
	yticklabelsvisible=false,
	yticksvisible=false,
	height=10,
	ygridvisible=false,
	xgridvisible=false,
)
linkxaxes!(metric_axes..., collisions_axis)

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
# Make coordinates and rotation responsive to the time slider
x, y, r = [@lift $data.robots[:, i, $(timestep.value)] for i in [X, Z, θ]]

# make color of robots dependent on clustering
c = @lift (
	if !$(robot_toggles[2].active)
		repeat([RGBA(0, 0, 0, 1)], size($data.robots, 1))
	else
		$discrete_palette[collect(
			cutree($data.clustering[$(timestep.value)]; h=$(threshold.value))
		)]
	end
)

# make glowcolor of robots dependend on collisions
g = @lift ( #TODO: refactor
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
				Makie.wong_colors()[6]
			elseif checkbounds(Bool, $wall_collisions, 1, $(timestep.value)) && any(
				$wall_collisions[
					i, max(($(timestep.value) - skip.value[]), 1):($(timestep.value))
				],
			)
				Makie.wong_colors()[7]
			else
				RGBA(0, 0, 0, 0)
			end for i in 1:size($wall_collisions, 1)
		]
	end
)

# Plot the robot swarm
#TODO: move center to center of mass
robot_marker = Makie.Polygon(Point2f[(-1, -1), (0, 0), (-1, 1), (2, 0)])
marker_scale = (@lift round(Int, 20 / sqrt(size($data.robots, 1))))
scatter!(
	swarm_animation,
	x,
	y;
	marker=robot_marker,
	markersize=marker_scale,
	rotations=r,
	color=c,
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
			if ndims(s) == 1
				lines!(axis, s; linewidth=1, color=Makie.wong_colors()[i], depth=1)
				limits!(
					axis,
					(nothing, nothing),
					maximum(s) <= 1 && minimum(s) >= 0 ? (0, 1) : (nothing, nothing),
				)
			else
				stride = 30
				transform = log
				for j in 1:(size(s, 1) + 1)
					band!(
						axis,
						(@lift 1:stride:($n_timesteps)),
						(@lift if j == 1
							transform.(
								fill(minimum(s[:, 1:stride:end]), $n_timesteps ÷ stride + 1)
							)
						else
							transform.(s[j, 1:stride:end])
						end),
						(@lift if j == (size(s, 1) + 1)
							transform.(s[j, 1:stride:end])
						else
							transform.(
								fill(maximum(s[:, 1:stride:end]), $n_timesteps ÷ stride + 1)
							)
						end);
						linewidth=1,
						depth_shift=-10000000,
					)
				end
				hlines!(
					axis,
					(@lift log($(threshold.value)));
					linewidth=0.5,
					color=:black,
					depth_shift=1000,
				)
			end
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
	) for (depth, (obs, color)) in enumerate(
		zip(
			[wall_collisions, agent_collisions],
			[Makie.wong_colors()[7], Makie.wong_colors()[6]],
		),
	)
]
for p in collision_plots
	connect!(p.visible, robot_toggles[1].active)
end

# plot timestep markers
for axis in metric_axes
	vlines!(axis, timestep.value; color=:black, linewidth=0.5, depth_shift=1000)
end
#TODO: add rectangle in collisions plot?

# Display the figure in it’s own window

display(fig);
