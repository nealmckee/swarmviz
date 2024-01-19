import AlgebraOfGraphics: set_aog_theme!
using Colors
using CSV
using DataFrames
using Distances
using GLMakie
using JSON
using LinearAlgebra
using NativeFileDialog
using NPZ
using Statistics

const ROBOTS = 1
const PROPERTIES = 2
const T = 3
const X = 2
const Y = 4
const θ = 5
const HVX = 6
const HVY = 7
const TRACKING_DIM = 5

include("src/metrics.jl")
include("src/plotstyle.jl")
include("src/analysis.jl")

"""
Struct that holds all the data input as well as everything we calculate.
"""
struct SwarmData
	"robots x properties x timesteps"
	tracking::Array{Float64,3}
	"robots x properties x timesteps"
	derived::Array{Float64,3}
	"metric name => datavector"
	analysis::Dict{String,Vector{Real}}
	"name => datavector"
	geometry::Dict{String,Any} # TODO: rename
end

# Set up observables
data = Observable(
	SwarmData(zeros(1, TRACKING_DIM, 1), zeros(1, TRACKING_DIM, 1), Dict(), Dict())
) #TODO: check empty dicts without dummy data
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(falses(1, 1))
agent_collisions = Observable(falses(1, 1))
n_timesteps = @lift size($data.tracking, T)
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

swarm_animation = Axis(fig[1, 1]; xlabel="X", ylabel="Y", autolimitaspect=1)

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
animation_toggles = [Toggle(fig) for _ in 1:4]
animation_controls[2:end, 3] = grid!(
	hcat(
		animation_toggles,
		[Label(fig, l; halign=:left) for l in ["Collisions", "Polygon", "CoM", "Diameter"]],
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
		l in ["Import Tracking", "Import Wall", "Import Collisions", "Export Metrics"]
	]

on(import_button.clicks) do c
	tracking_path = pick_file(; filterlist="npy")
	tracking_path != "" && (data[] = analyse_tracking(tracking_path))
	autolimits!(metric_axes[1])
	limits!(
		swarm_animation,
		minimum(data[].tracking[:, 2, :]) - 100,
		maximum(data[].tracking[:, 2, :]) + 100,
		minimum(data[].tracking[:, 4, :]) - 100,
		maximum(data[].tracking[:, 4, :]) + 100,
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
	CSV.write(joinpath(export_folder, "metrics.csv"), DataFrame(data[].analysis))
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
metrics_tuples = @lift Tuple.(collect($data.analysis))
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
x, y, r = [@lift $data.tracking[:, i, $(time_slider.sliders[1].value)] for i in [2, 4, 5]]

# make color of robots dependend on collisions
g = @lift ( #TODO: refactor
	if size($agent_collisions, 1) != size($data.tracking, 1) ||
		size($wall_collisions, 1) != size($data.tracking, 1) ||
		!$(animation_toggles[1].active)
		repeat([RGBA(0, 0, 0, 0)], size($data.tracking, 1))
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
	color=:black,
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
	@lift Point2f.($data.geometry["Surrounding Polygon"][$(time_slider.sliders[1].value)]);
	color=:transparent,
	strokecolor="#8f8f8f",
	strokewidth=1,
	linestyle=:dash,
	closed=true,
)
connect!(surrounding_polygon.visible, animation_toggles[2].active)

# plot the center of mass and connect to toggle
center_of_mass = scatter!(
	swarm_animation,
	@lift Point2f($data.geometry["Center of Mass"][$(time_slider.sliders[1].value)]);
	color="#8f8f8f",
	markersize=12,
	marker=:xcross,
)
connect!(center_of_mass.visible, animation_toggles[3].active)

#plot the diameter and connect to toggle
diameter = lines!(
	swarm_animation,
	@lift Point2f.(
		eachslice(
			$data.tracking[
				[
					Tuple(
						$data.geometry["Furthest Robots"][$(time_slider.sliders[1].value)]
					)...,
				],
				[X, Y],
				$(time_slider.sliders[1].value),
			];
			dims=ROBOTS,
		)
	);
	color="#8f8f8f",
	linewidth=1,
)
connect!(diameter.visible, animation_toggles[4].active)

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
			limits!(axis, (nothing, nothing), maximum(s) <= 1 ? (0, 1) : (nothing, nothing))
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
	connect!(p.visible, animation_toggles[1].active)
end

# plot timestep markers
for axis in metric_axes
	vlines!(axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
end
#TODO: add rectangle in collisions plot?

# Display the figure in it’s own window
fig
