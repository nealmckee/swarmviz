using AlgebraOfGraphics
using Colors
using Distances
using GLMakie
using LinearAlgebra
using NativeFileDialog
using NPZ
using Statistics

wall_data = zeros(2, 1)
tracking_data = zeros(1, 7, 1)
metrics_data = zeros(3, 1)

include("src/metrics.jl")
include("src/plotstyle.jl")

experiment_file = "data/DatasetE2/E21/E212/E212r1_summaryd.npy"
wall_file = "data/DatasetE2/ArenaBorders/ArenaBorders_r1_summaryd.npy"

# Load the data and drop singular dimensions
tracking_data = npzread(experiment_file)[1, :, :, :]
wall_data = npzread(wall_file)[1, 1, [2, 4], :]

# Remove columns with NaN values from wall_data
wall_data = wall_data[:, .!isnan.(wall_data[1, :])]
wall_data = wall_data[:, .!isnan.(wall_data[2, :])]

# Get the number of robots and the number of timesteps
n_robots = size(tracking_data, 1)
n_timesteps = size(tracking_data, 3)

# Add heading vector from orientation to data
tracking_data = cat(
    tracking_data, reshape(cos.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps); dims=2
)
tracking_data = cat(
    tracking_data, reshape(sin.(tracking_data[:, 5, :]), n_robots, 1, n_timesteps); dims=2
)

# precalculate all metrics for all timesteps (currently no clustering)
polarisation = dropdims(
    mapslices(swarm_polarisation, tracking_data; dims=(1, 2)); dims=(1, 2)
)
rotational_order = dropdims(
    mapslices(swarm_rotational_order, tracking_data; dims=(1, 2)); dims=(1, 2)
)
mean_interindividual_distance = dropdims(
    mapslices(swarm_mean_interindividual_distance, tracking_data; dims=(1, 2)); dims=(1, 2)
)
metrics_data = cat(polarisation', rotational_order', mean_interindividual_distance'; dims=1)

# Set up the figure
GLMakie.activate!(; title="SwarmViz")
fig = Figure(; size=(960, 600))

swarm_animation = Axis(fig[1:3, 1:2]; xlabel="X", ylabel="Y")

time_slider = SliderGrid(
    fig[4, 1:2], (label="Timestep", range=1:1:n_timesteps, startvalue=1)
)

video_control = Button(fig[5, 1]; label="Play/Pause")

video_settings = SliderGrid(
    fig[5, 2],
    (label="FPS", range=1:1:240, startvalue=30, format="{:.1f}Hz"),
    (label="Skip", range=0:1:240, startvalue=0),
)

isplaying = Observable(false)
on(video_control.clicks) do c
    isplaying[] = !isplaying[]
end
on(video_control.clicks) do c
    @async while isplaying[] && #TODO: check whether @async is safe/necessary
                 time_slider.sliders[1].value[] <
                 n_timesteps - video_settings.sliders[2].value[] - 1
        set_close_to!(
            time_slider.sliders[1],
            time_slider.sliders[1].value[] + 1 + video_settings.sliders[2].value[],
        )
        sleep(1 / video_settings.sliders[1].value[])
        isopen(fig.scene) || break # crucial, ensures computations stop if closed window.
    end
end

polarisation_axis = Axis(
    fig[1, 3:4]; ylabel="Polarisation", limits=(nothing, (0, 1)), xticklabelsvisible=false
)
rotational_order_axis = Axis(
    fig[2, 3:4];
    ylabel="Rotational Order",
    limits=(nothing, (0, 1)),
    xticklabelsvisible=false,
)
mean_interindividual_distance_axis = Axis(fig[3, 3:4]; ylabel="Mean IID", xlabel="Timestep")

buttongrid = GridLayout(fig[4:5, 3]; default_rowgap=4)

import_button, wall_button, export_button =
    buttongrid[1:3, 1] = [
        Button(fig; label=l, halign = :left) for l in ["Import Tracking", "Import Wall", "Export Analysis"]
    ]


on(import_button.clicks) do c
    experiment_file = pick_file(; filterlist="*.npy")
end

on(wall_button.clicks) do c
    wall_file = pick_file(; filterlist="*.npy")
end

# Plot the wall of the enclosure
poly!(
    swarm_animation,
    Point2f.(wall_data[1, :], wall_data[2, :]);
    color=:transparent,
    strokecolor="#8f8f8f",
    strokewidth=1,
    linestyle=:dot,
    closed=true,
)

# Make coordinates and heading vectors responsive to the time slider
x, y, u, v = [
    lift(time_slider.sliders[1].value) do t
        tracking_data[:, i, t]
    end for i in [2, 4, 6, 7]
]

# Plot the robot swarm
arrows!(swarm_animation, x, y, u, v; lengthscale=100)
scatter!(swarm_animation, x, y; markersize=12, color=:black)

# Plot the metrics #TODO: for loop
lines!(polarisation_axis, 1:n_timesteps, metrics_data[1, :]; linewidth=1)
vlines!(polarisation_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(rotational_order_axis, 1:n_timesteps, metrics_data[2, :]; linewidth=1)
vlines!(rotational_order_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(mean_interindividual_distance_axis, 1:n_timesteps, metrics_data[3, :]; linewidth=1)
vlines!(
    mean_interindividual_distance_axis,
    time_slider.sliders[1].value;
    color=:black,
    linewidth=0.5,
)

# Display the figure in it’s own window
fig
