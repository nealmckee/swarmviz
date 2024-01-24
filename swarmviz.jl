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
	metrics::Dict{String,Vector{Float64}}
	"name => vector with derived ata for each timestep"
	derived::Dict{String,Any}
	"one clustering per timesteps"
	clustering::Vector{Hclust{Float64}}
end

# Set up observables
#TODO: check empty dicts without dummy data
data = Observable(SwarmData(zeros(1, TRACKING_DIM + 10, 1), Dict(), Dict(), []))
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(falses(1, 1))
agent_collisions = Observable(falses(1, 1))
n_timesteps = @lift size($data.robots, T)
timesteps = @lift 1:($n_timesteps)
isplaying = Observable(false)
time()

# Preload data for easier debugging #TODO: remove when done
tracking_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_summaryd.npy"
wall_path = "data/A0_09_B0_0/w3_summaryd.npy"
wall_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_warefl.json"
agent_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_aarefl.json"
wall_data[] = analyse_wall(wall_path)
data[] = analyse_tracking(tracking_path)
wall_collisions[] = process_collisions(wall_collisons_path)
agent_collisions[] = process_collisions(agent_collisons_path)

# Set up the figure
GLMakie.activate!(; title="SwarmViz")
fig = Figure(; size=(960, 600))

swarm_animation = Axis(fig[1, 1]; xlabel="X", ylabel="Z", autolimitaspect=1)

#TODO: use space on the left
time_slider = SliderGrid(fig[2, 1:2], (label="Timestep", range=timesteps, startvalue=1))
animation_controls = GridLayout(fig[3, 1]; default_rowgap=12, default_colgap=12)

metrics_grid = GridLayout(fig[1, 2])
data_controls = GridLayout(fig[3, 2]; default_rowgap=2)

#adjust the width of the column with the swarm animation
colsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))

# create play/pause button with reactive label
play_button_text = @lift $isplaying[] ? "Pause" : "Play"
play_button = Button(animation_controls[2, 1]; label=play_button_text, width=60)

animation_settings = SliderGrid(
	animation_controls[2, 2],
	(label="FPS", range=1:1:120, startvalue=30),
	(label="Skip", range=0:1:99, startvalue=0),
)
animation_toggles = [Toggle(fig) for _ in 1:5]
animation_controls[2:end, 3] = grid!(
	hcat(
		animation_toggles,
		[
			Label(fig, l; halign=:left) for
			l in ["Clustering", "Collisions", "Polygon", "CoM", "Diameter"]
		],
	);
	default_rowgap=3,
	default_colgap=3,
)

# start animation loop on buttonpress, the plot then automatically updates
# as it’s dependent on the time slider
on(play_button.clicks) do c
	isplaying[] = !isplaying[]
end
on(play_button.clicks) do c
	@async while isplaying[] &&
				 time_slider.sliders[1].value[] <
				 n_timesteps[] - animation_settings.sliders[2].value[] - 1
		frame_start = time()
		set_close_to!(
			time_slider.sliders[1],
			time_slider.sliders[1].value[] + 1 + animation_settings.sliders[2].value[],
		)
		sleep(max(1 / animation_settings.sliders[1].value[] - time() + frame_start, 0))
		isopen(fig.scene) || break # ensures computations stop if window is closed
		if time_slider.sliders[1].value[] >=
			n_timesteps[] - animation_settings.sliders[2].value[] - 1
			isplaying[] = !isplaying[]
		end
	end
end

# create buttons to import/export data
import_button, wall_button, collision_button, export_metrics_button =
	data_controls[1:2, 1:2] = [
		Button(fig; label=l, halign=:left) for
		l in ["Import Tracking", "Import Wall", "Import Collisions", "Export All"]
	]

on(import_button.clicks) do c
	tracking_path = pick_file(; filterlist="npy")
	tracking_path != "" && (data[] = analyse_tracking(tracking_path))
	autolimits!(metric_axes[1])
	limits!(
		swarm_animation,
		minimum(data[].robots[:, X, :]) - 100,
		maximum(data[].robots[:, X, :]) + 100,
		minimum(data[].robots[:, Z, :]) - 100,
		maximum(data[].robots[:, Z, :]) + 100,
	)
	set_close_to!(time_slider.sliders[1], 1)
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
metric_axes = [Axis(metrics_grid[row, 1]; xticklabelsvisible=false) for row in 1:3]
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
metrics_tuples = @lift Tuple.(collect($data.metrics))
metric_menus = [
	Menu(
		metrics_grid[i, 1, Top()];
		options=metrics_tuples,
		default=default,
		width=150,
		height=18,
		halign=:center,
		prompt="Select Metric...",
	) for (i, default) in enumerate(["Polarisation", "Rotational Order", "Diameter"]) #TODO remove defaults after debugging
]
# Make coordinates and rotation responsive to the time slider
x, y, r = [@lift $data.robots[:, i, $(time_slider.sliders[1].value)] for i in [X, Z, θ]]

# make color of robots dependent on clustering
c = @lift (
	if !$(animation_toggles[1].active)
		repeat([RGBA(0, 0, 0, 1)], size($data.robots, 1))
	else
		Makie.wong_colors()[collect(
			cutree($data.clustering[$(time_slider.sliders[1].value)]; k=3)
		)]
	end
)

# make glowcolor of robots dependend on collisions
g = @lift ( #TODO: refactor
	if size($agent_collisions, 1) != size($data.robots, 1) ||
		size($wall_collisions, 1) != size($data.robots, 1) ||
		!$(animation_toggles[2].active)
		repeat([RGBA(0, 0, 0, 0)], size($data.robots, 1))
	else
		[
			if checkbounds(Bool, $agent_collisions, 1, $(time_slider.sliders[1].value)) &&
				any(
				$agent_collisions[
					i,
					max(
						(
							$(time_slider.sliders[1].value) -
							animation_settings.sliders[2].value[]
						),
						1,
					):($(time_slider.sliders[1].value)),
				],
			)
				Makie.wong_colors()[6]
			elseif checkbounds(
				Bool, $wall_collisions, 1, $(time_slider.sliders[1].value)
			) && any(
				$wall_collisions[
					i,
					max(
						(
							$(time_slider.sliders[1].value) -
							animation_settings.sliders[2].value[]
						),
						1,
					):($(time_slider.sliders[1].value)),
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
scatter!(
	swarm_animation,
	x,
	y;
	marker=robot_marker,
	markersize=6,
	rotations=r,
	color=c,
	glowcolor=g,
	glowwidth=6,
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
	@lift Point2f.($data.derived["Convex Hull"][$(time_slider.sliders[1].value)]);
	color=:transparent,
	strokecolor="#8f8f8f",
	strokewidth=1,
	linestyle=:dash,
	closed=true,
)
connect!(surrounding_polygon.visible, animation_toggles[3].active)

# plot the center of mass and connect to toggle
center_of_mass = scatter!(
	swarm_animation,
	@lift Point2f($data.derived["Center of Mass"][$(time_slider.sliders[1].value)]);
	color="#8f8f8f",
	markersize=12,
	marker=:xcross,
)
connect!(center_of_mass.visible, animation_toggles[4].active)

#plot the diameter and connect to toggle
diameter = lines!(
	swarm_animation,
	@lift Point2f.(
		eachslice(
			$data.robots[
				$data.derived["Furthest Robots"][$(time_slider.sliders[1].value)],
				[X, Z],
				$(time_slider.sliders[1].value),
			];
			dims=ROBOTS,
		)
	);
	color="#8f8f8f",
	linewidth=1,
)
connect!(diameter.visible, animation_toggles[5].active)

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
			lines!(axis, s; linewidth=1, color=Makie.wong_colors()[i])
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
	) for (depth, (obs, color)) in enumerate(
		zip(
			[wall_collisions, agent_collisions],
			[Makie.wong_colors()[7], Makie.wong_colors()[6]],
		),
	)
]
for p in collision_plots
	connect!(p.visible, animation_toggles[2].active)
end

# plot timestep markers
for axis in metric_axes
	vlines!(axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
end
#TODO: add rectangle in collisions plot?

# Display the figure in it’s own window
display(fig)
