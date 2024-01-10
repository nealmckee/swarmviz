import AlgebraOfGraphics: set_aog_theme!
using Colors
using Distances
using GLMakie
using JSON
using LinearAlgebra
using NativeFileDialog
using NPZ
using Statistics

include("src/metrics.jl")
include("src/plotstyle.jl")
include("src/analysis.jl")

struct SwarmData
	tracking::Array{Float64,3} # robots x properties x timesteps
	derived::Array{Float64,3} # robots x properties x timesteps
	analysis::Dict{String,Any} # name, datamatrix
	geometry::Dict{String,Any} # name, datavector
end

# Set up observables
data = Observable(SwarmData(zeros(1, 5, 1), zeros(1, 5, 1), Dict(), Dict()))
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(falses(1, 1))
agent_collisions = Observable(falses(1, 1))
n_timesteps = @lift size($data.tracking, 3)
timesteps = @lift 1:($n_timesteps)

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
animation_controls = GridLayout(fig[2, 1]; default_rowgap=12, default_colgap=12)
time_slider = SliderGrid(
	animation_controls[1, 1:3], (label="Timestep", range=timesteps, startvalue=1)
)

# Watch for Play/Pause status
isplaying = Observable(false)
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

# Watch for changes in the time slider to update the plot
on(play_button.clicks) do c
	isplaying[] = !isplaying[]
end
on(play_button.clicks) do c
	@async while isplaying[] && #TODO: check whether @async is safe/necessary
				 time_slider.sliders[1].value[] <
				 n_timesteps[] - animation_settings.sliders[2].value[] - 1
		set_close_to!(
			time_slider.sliders[1],
			time_slider.sliders[1].value[] + 1 + animation_settings.sliders[2].value[],
		)
		sleep(1 / animation_settings.sliders[1].value[])
		isopen(fig.scene) || break # ensures computations stop if window is closed
	end
end

metrics_grid = GridLayout(fig[1, 2])
metric_axis1 = Axis(metrics_grid[1, 1]; xticklabelsvisible=false)
metric_axis2 = Axis(metrics_grid[2, 1]; xticklabelsvisible=false)
metric_axis3 = Axis(metrics_grid[3, 1]; xticklabelsvisible=false)
collisions_axis = Axis(
	metrics_grid[4, 1];
	xlabel="Timestep",
	yticklabelsvisible=false,
    yticksvisible=false,
	height=10,
	ygridvisible=false,
    xgridvisible=false,
)
linkxaxes!(metric_axis1, metric_axis2, metric_axis3, collisions_axis)

data_controls = GridLayout(fig[2, 2]; default_rowgap=2)
buttongrid = GridLayout(data_controls[1:3, :]; default_rowgap=2)

import_button, wall_button, collision_button =
	buttongrid[1:3, 1] = [
		Button(fig; label=l, halign=:left) for
		l in ["Import Tracking", "Import Wall", "Import Collisions"]
	]

on(import_button.clicks) do c
	tracking_path = pick_file(; filterlist="npy")
	tracking_path != "" && (data[] = analyse_tracking(tracking_path))
	autolimits!(metric_axis3)
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
end

# menus to select which metrics to plot
metric_menu1 = Menu(
	metrics_grid[1, 1, Top()];
	options=Tuple.(collect(data[].analysis)),
	default="Polarisation",
	width=150,
	height=18,
	halign=:center,
	prompt="Select Metric...",
)
metric_menu2 = Menu(
	metrics_grid[2, 1, Top()];
	options=Tuple.(collect(data[].analysis)),
	default="Rotational Order",
	width=150,
	height=18,
	halign=:center,
	prompt="Select Metric...",
)
metric_menu3 = Menu(
	metrics_grid[3, 1, Top()];
	options=Tuple.(collect(data[].analysis)),
	default="Mean IID",
	width=150,
	height=18,
	halign=:center,
	prompt="Select Metric...",
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

# Make coordinates and rotation responsive to the time slider
x, y, r = [@lift $data.tracking[:, i, $(time_slider.sliders[1].value)] for i in [2, 4, 5]]

# make color of robots dependend on collisions
c = @lift ( #TODO: refactor
	if size($agent_collisions, 1) != size($data.tracking, 1) ||
		size($wall_collisions, 1) != size($data.tracking, 1) ||
		!$(animation_toggles[1].active)
		repeat([RGBA(0, 0, 0)], size($data.tracking, 1))
	else
		[
			if checkbounds(Bool, $agent_collisions, 1, $(time_slider.sliders[1].value)) &&
				any(
				$agent_collisions[
					i,
					($(time_slider.sliders[1].value) - animation_settings.sliders[2].value[]):($(
						time_slider.sliders[1].value
					)),
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
				RGBA(0, 0, 0)
			end for i in 1:size($wall_collisions, 1)
		]
	end
)

# Plot the robot swarm
robot_marker = Makie.Polygon(Point2f[(-1, -1), (0, 0), (-1, 1), (2, 0)])
#TODO: move center to center of mass
#TODO: switch collisions to glow instead of color?
scatter!(swarm_animation, x, y; marker=robot_marker, markersize=6, rotations=r, color=c)

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
center_of_mass = scatter!( #TODO: plot styling (feedback?)
	swarm_animation,
	@lift Point2f($data.geometry["Center of Mass"][$(time_slider.sliders[1].value)]);
	color="#8f8f8f",
	markersize=12,
	marker=:xcross,
)
connect!(center_of_mass.visible, animation_toggles[3].active)

#plot the diameter and connect to toggle
diameter = lines!( #TODO: plot styling (feedback?)
	swarm_animation,
	@lift Point2f.(
		$data.geometry["Surrounding Polygon"][$(time_slider.sliders[1].value)][[
			Tuple($data.geometry["Furthest Robots"][$(time_slider.sliders[1].value)])...
		]]
	);
	color="#8f8f8f",
	linewidth=1,
)
connect!(diameter.visible, animation_toggles[4].active)

# plot the metrics when one is chosen from the corresponding menu
for (i, (menu, axis)) in enumerate(
	zip(
		[metric_menu1, metric_menu2, metric_menu3],
		[metric_axis1, metric_axis2, metric_axis3],
	),
)
	on(menu.selection) do s
		while length(axis.scene.plots) > 1
			delete!(
				axis,
				axis.scene.plots[typeof.(axis.scene.plots) .!= Plot{
					Makie.vlines,Tuple{Int64}
				}][1],
			)
		end
		lines!(axis, s; linewidth=1, color=Makie.wong_colors()[i])
		limits!(axis, (nothing, nothing), maximum(s) <= 1 ? (0, 1) : (nothing, nothing))
	end
	notify(menu.selection)
end

# plot timestep markers
for axis in [metric_axis1, metric_axis2, metric_axis3]
	vlines!(axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
end
#TODO: add rectangle in collisions plot?

# plot collisions
agent_collisions_plot = vlines!(
    collisions_axis,
    @lift findall(
        reduce(|, $agent_collisions, dims = 1)[1,:]
    );
    color=Makie.wong_colors()[6],
    linewidth=0.1,
    depth = 2,
)
wall_collisions_plot = vlines!(
    collisions_axis,
    @lift findall(
        reduce(|, $wall_collisions, dims = 1)[1,:]
    );
    color=Makie.wong_colors()[7],
    linewidth=0.1,
    depth = 1,
)
connect!(agent_collisions_plot.visible, animation_toggles[1].active)
connect!(wall_collisions_plot.visible, animation_toggles[1].active)

#adjust the width of the column with the swarm animation
colsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))

# Display the figure in it’s own window
fig
