using AlgebraOfGraphics
using Colors
using Distances
using GLMakie
using LinearAlgebra
using NPZ
using Statistics

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

# Calculate the polarisation of the swarm for the data of one timestep
function swarm_polarisation(tracking_datapoint; cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    return norm(sum(tracking_datapoint[:, 6:7]; dims=1)) / size(tracking_datapoint, 1)
end

# Calculate the rotational order of the swarm for the data of one timestep
function swarm_rotational_order(tracking_datapoint; cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    clustersize = size(tracking_datapoint, 1)
    barycenter = sum(tracking_datapoint[:, [2, 4]]; dims=1) ./ clustersize
    radial_unit_vectors = mapslices(
        x -> x / norm(x), tracking_datapoint[:, [2, 4]] .- barycenter; dims=2
    )
    return norm(
        sum([
            cross(
                push!(tracking_datapoint[i, 6:7], 0), push!(radial_unit_vectors[i, :], 0)
            ) for i in 1:clustersize
        ]),
    ) / clustersize
end

function swarm_mean_interindividual_distance(tracking_datapoint, cluster_members=[])
    if !isempty(cluster_members)
        tracking_datapoint = tracking_datapoint[cluster_members, :]
    end
    return mean(pairwise(Euclidean(), tracking_datapoint[:, [2, 4]]; dims=1))
end

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

metrics_data = cat(polarisation, rotational_order, mean_interindividual_distance; dims=2)

# Plot styling
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
        # xticklabelfont=lightfont,
        # yticklabelfont=lightfont,
        # xlabelfont=mediumfont,
        # ylabelfont=mediumfont,
        # titlefont=mediumfont,
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
        # protrusions=30, # seems to be fixed in Makie
        #     # xticklabelfont=lightfont,
        #     # yticklabelfont=lightfont,
        #     # zticklabelfont=lightfont,
        #     # xlabelfont=mediumfont,
        #     # ylabelfont=mediumfont,
        #     # zlabelfont=mediumfont,
        #     # titlefont=mediumfont,
    ),
    Legend=(
        framevisible=false,
        gridshalign=:left,
        padding=(0.0f0, 0.0f0, 0.0f0, 0.0f0),
        labelcolor="#8f8f8f",
        titlecolor="#8f8f8f",
        # labelfont=lightfont,
        # titlefont=mediumfont,
    ),
    Colorbar=(
        tickcolor="#8f8f8f",
        #     flip_vertical_label = true,
        #     spinewidth = 0,
        #     # ticklabelfont=lightfont,
        #     # labelfont=mediumfont,
    ),
    resolution=(800, 500),
)

# Set up the figure
GLMakie.activate!(; title="SwarmViz")
fig = Figure(; resolution=(960, 600))

swarm_animation = Axis(fig[1:3, 1:2]; xlabel="X", ylabel="Y")

time_slider = SliderGrid(
    fig[4, 1:2], (label="Timestep", range=1:1:n_timesteps, startvalue=1)
)

video_settings = SliderGrid(
    fig[5, 2],
    (label="FPS", range=1:1:240, startvalue=30, format="{:.1f}Hz"),
    (label="Skip", range=0:1:240, startvalue=0),
)

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
x = lift(time_slider.sliders[1].value) do timestep
    tracking_data[:, 2, timestep]
end
y = lift(time_slider.sliders[1].value) do timestep
    tracking_data[:, 4, timestep]
end
u = lift(time_slider.sliders[1].value) do timestep
    tracking_data[:, 6, timestep]
end
v = lift(time_slider.sliders[1].value) do timestep
    tracking_data[:, 7, timestep]
end

# Plot the robot swarm
arrows!(swarm_animation, x, y, u, v; lengthscale=100)

# Plot the metrics
lines!(polarisation_axis, 1:n_timesteps, polarisation; linewidth=1)
vlines!(polarisation_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(rotational_order_axis, 1:n_timesteps, rotational_order; linewidth=1)
vlines!(rotational_order_axis, time_slider.sliders[1].value; color=:black, linewidth=0.5)
lines!(
    mean_interindividual_distance_axis,
    1:n_timesteps,
    mean_interindividual_distance;
    linewidth=1,
)
vlines!(
    mean_interindividual_distance_axis,
    time_slider.sliders[1].value;
    color=:black,
    linewidth=0.5,
)

# Display the figure in it’s own window
fig
