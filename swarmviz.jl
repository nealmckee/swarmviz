using AlgebraOfGraphics
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
	analysis::Array{Float64,2} # metrics x timesteps #TODO: better typing? another struct?
end

# Set up observables
data = Observable(SwarmData(zeros(1, 5, 1), zeros(1, 5, 1), zeros(3, 1)))
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(BitArray(undef, 1, 1))
agent_collisions = Observable(BitArray(undef, 1, 1))
n_timesteps = @lift size($data.tracking, 3)
timesteps = @lift 1:($n_timesteps)

# Preload data for easier debugging #TODO: remove when done
tracking_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_summaryd.npy"
wall_path = "data/DatasetE2/ArenaBorders/ArenaBorders_r1_summaryd.npy"
wall_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_warefl.json"
agent_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_aarefl.json"
wall_data[] = analyse_wall(wall_path)
data[] = analyse_tracking(tracking_path)
wall_collisions[] = process_collisions(wall_collisons_path)
agent_collisions[] = process_collisions(agent_collisons_path)

# Set up the figure #TODO make layout more organised and fix scaling of blocks
GLMakie.activate!(; title="SwarmViz")
fig = Figure(; size=(960, 600))

swarm_animation = Axis(fig[1, 1]; xlabel="X", ylabel="Y", autolimitaspect=1) #TODO: fix aspect ratio
animation_controls = GridLayout(fig[2, 1]; tellwidth=false)
time_slider = SliderGrid(
	animation_controls[1, 1:2], (label="Timestep", range=timesteps, startvalue=1);
)

# Watch for Play/Pause status
isplaying = Observable(false)
play_button_text = @lift $isplaying[] ? "Pause" : "Play"
play_button = Button(animation_controls[2, 1]; label=play_button_text, width=60)

animation_settings = SliderGrid(
	animation_controls[2, 2],
	(label="FPS", range=1:1:240, startvalue=30, format="{:.1f}Hz"),
	(label="Skip", range=0:1:240, startvalue=0);
	tellwidth=false,
)

# Watch for changes in the time slider to update the plot
on(play_button.clicks) do c
	return isplaying[] = !isplaying[]
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
pol_axis = Axis(
	metrics_grid[1, 1];
	ylabel="Polarisation",
	limits=(nothing, (0, 1)),
	xticklabelsvisible=false,
)
ro_axis = Axis(
	metrics_grid[2, 1];
	ylabel="Rotational Order",
	limits=(nothing, (0, 1)),
	xticklabelsvisible=false,
)
miid_axis = Axis(metrics_grid[3, 1]; ylabel="Mean IID", xlabel="Timestep")
linkxaxes!(pol_axis, ro_axis, miid_axis)

file_controls = GridLayout(fig[2, 2])
buttongrid = GridLayout(file_controls[1, 1]; default_rowgap=4)

import_button, wall_button, export_button =
	buttongrid[1:3, 1] = [
		Button(fig; label=l, halign=:left) for
		l in ["Import Tracking", "Import Wall", "Export Analysis"]
	]

on(import_button.clicks) do c
	tracking_path = pick_file(; filterlist="npy")
	tracking_path != "" && (data[] = analyse_tracking(tracking_path))
	autolimits!(miid_axis)
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
)

# Make coordinates and rotation responsive to the time slider
x, y, r = [@lift $data.tracking[:, i, $(time_slider.sliders[1].value)] for i in [2, 4, 5]]

# Plot the robot swarm
robot_marker = Makie.Polygon(Point2f[(-1, -1), (0, 0), (-1, 1), (2, 0)])
scatter!(
	swarm_animation,
	x,
	y;
	marker=robot_marker,
	markersize=6,
	rotations=r,
	color=(@lift [
        if $agent_collisions[i, $(time_slider.sliders[1].value)]
            Makie.wong_colors()[6]
		elseif $wall_collisions[i, $(time_slider.sliders[1].value)]
			Makie.wong_colors()[7]
		else
			RGBA(0,0,0)
		end for i in 1:size($wall_collisions, 1)
	]),
)

# Plot the metrics #TODO: for loop
lines!(
	pol_axis, @lift float.($data.analysis[1, :]); linewidth=1, color=Makie.wong_colors()[1]
)
vlines!(pol_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(
	ro_axis, @lift float.($data.analysis[2, :]); linewidth=1, color=Makie.wong_colors()[2]
)
vlines!(ro_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(
	miid_axis, @lift float.($data.analysis[3, :]); linewidth=1, color=Makie.wong_colors()[3]
)
vlines!(miid_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)

colsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))

# Display the figure in it’s own window
fig
