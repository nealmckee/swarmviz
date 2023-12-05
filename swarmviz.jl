using AlgebraOfGraphics
using Colors
using Distances
using GLMakie
using LinearAlgebra
using NativeFileDialog
using NPZ
using Statistics

include("src/metrics.jl")
include("src/plotstyle.jl")
include("src/analysis.jl")

struct SwarmData
    tracking::Array{Float64,3} # robots x datapoints x timesteps
    derived::Array{Float64,3}
    analysis::Array{Float64,2} # metrics x timesteps #TODO: better typing? another struct?
end

# Set up observables
data = Observable(SwarmData(zeros(1, 5, 1), zeros(1, 5, 1), zeros(3, 1)))
wall_data = Observable(zeros(2, 1))
n_timesteps = @lift size($data.tracking, 3)
timesteps = @lift 1:($n_timesteps)

# Preload data for easier debugging #TODO: remove when done
tracking_path = "data/DatasetE2/E22/E223/E223r1_summaryd.npy"
wall_path = "data/DatasetE2/ArenaBorders/ArenaBorders_r1_summaryd.npy"
wall_data[] = analyse_wall(wall_path)
data[] = analyse_tracking(tracking_path)

# Set up the figure #TODO make layout more organised and fix scaling of blocks
GLMakie.activate!(; title="SwarmViz")
fig = Figure(; size=(960, 600))

swarm_animation = Axis(fig[1:3, 1:2]; xlabel="X", ylabel="Y", aspect=1) #TODO: fix aspect ratio

time_slider = SliderGrid(
    fig[4, 1:2], (label="Timestep", range=timesteps, startvalue=1); tellwidth=false
)

# Watch for Play/Pause status
isplaying = Observable(false)
play_button_text = @lift $isplaying[] ? "Pause" : "Play"
video_control = Button(fig[5, 1]; label=play_button_text, width=60)

video_settings = SliderGrid(
    fig[5, 2],
    (label="FPS", range=1:1:240, startvalue=30, format="{:.1f}Hz"),
    (label="Skip", range=0:1:240, startvalue=0),
)

# Watch for changes in the time slider to update the plot
on(video_control.clicks) do c
    return isplaying[] = !isplaying[]
end
on(video_control.clicks) do c
    @async while isplaying[] && #TODO: check whether @async is safe/necessary
                 time_slider.sliders[1].value[] <
                 n_timesteps[] - video_settings.sliders[2].value[] - 1
        set_close_to!(
            time_slider.sliders[1],
            time_slider.sliders[1].value[] + 1 + video_settings.sliders[2].value[],
        )
        sleep(1 / video_settings.sliders[1].value[])
        isopen(fig.scene) || break # ensures computations stop if window is closed
    end
end

pol_axis = Axis(
    fig[1, 3:4]; ylabel="Polarisation", limits=(nothing, (0, 1)), xticklabelsvisible=false
)
ro_axis = Axis(
    fig[2, 3:4];
    ylabel="Rotational Order",
    limits=(nothing, (0, 1)),
    xticklabelsvisible=false,
)
miid_axis = Axis(fig[3, 3:4]; ylabel="Mean IID", xlabel="Timestep")

buttongrid = GridLayout(fig[4:5, 3]; default_rowgap=4)

import_button, wall_button, export_button =
    buttongrid[1:3, 1] = [
        Button(fig; label=l, halign=:left) for
        l in ["Import Tracking", "Import Wall", "Export Analysis"]
    ]

on(import_button.clicks) do c
    tracking_path = pick_file(; filterlist="npy")
    tracking_path != "" && (data[] = analyse_tracking(tracking_path))
    autolimits!.([swarm_animation, pol_axis, ro_axis, miid_axis])
    set_close_to!(time_slider.sliders[1], 1)
end

on(wall_button.clicks) do c
    wall_path = pick_file(; filterlist="npy")
    wall_path == "" || (wall_data = analyse_wall(wall_path))
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
    swarm_animation, x, y; marker=robot_marker, markersize=6, rotations=r, color=:black
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

# Display the figure in it’s own window
fig
