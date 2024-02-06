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

# constants for indexing as specified for the input data format in the readme
const ROBOTS = 1
const PROPERTIES = 2
const TIME = 3
const T = 1
const X = 2
const Z = 4
const θ = 5
# TODO: add more constants if we require denoised derivatives
const HVX = 6
const HVZ = 7
const ACCX = 11 #TODO remove?
const ACCZ = 12 #TODO remove?
const TRACKING_DIM = 5 #TODO change if more inputs required or make adaptive?
# operating system dependent path to the assets folder
FONTFOLDER = RelocatableFolders.@path joinpath(@__DIR__, "assets") #TODO make const
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
include("src/metrics.jl")
include("src/analysis.jl")
include("src/dataprocessing.jl")
# include("src/tvdiff.jl") #TODO remove after final rejection

"""
Struct that holds all the data input as well as everything we calculate.
"""
struct SwarmData
	"robots x properties x timesteps"
	robots::Array{Float64,3}
	"metric name => vector with metric values for each timestep"
	metrics::Dict{String,Vector{Float64}}
	"name => vector with derived data for each timestep"
	derived::Dict{String,Any}
	"one clustering object per timesteps"
	clustering::Vector{Hclust{Float64}}
end

# Set up dummy observables
#TODO: check without dummy data
#TODO: Clustering can’t be empty
data = Observable(SwarmData(zeros(1, TRACKING_DIM + 10, 1), Dict(), Dict(), []))
wall_data = Observable(zeros(2, 1))
wall_collisions = Observable(falses(1, 1))
agent_collisions = Observable(falses(1, 1))
n_timesteps = @lift size($data.robots, TIME)
timesteps = @lift 1:($n_timesteps)
isplaying = Observable(false)
discrete_palette = Observable(PALETTE)

# Preload dummy data for easier debugging #TODO: remove when done
tracking_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_summaryd.npy"
wall_path = "data/A0_09_B0_0/w3_summaryd.npy"
wall_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_warefl.json"
agent_collisons_path = "data/A0_09_B0_0/EXP1_A0_09_B0_0_r1_w3_aarefl.json"
wall_data[] = analyse_wall(wall_path)
data[] = analyse_tracking(tracking_path)
wall_collisions[] = process_collisions(wall_collisons_path)
agent_collisions[] = process_collisions(agent_collisons_path)
discrete_palette[] = vcat(
	RGB.(PALETTE),
	distinguishable_colors(
		size(data[].robots, 1) - length(PALETTE),
		vcat(RGB.(PALETTE), [RGB(0, 0, 0), RGB(1, 1, 1)]);
		dropseed=true,
		lchoices=range(45, 90; length=15),
	),
)

# run the scripts to set up the style and layout,
# add interactivity, and create reactive plots
#! Order of execution matters and should never be changed!
include("scripts/plotstyle.jl")
include("scripts/layout.jl")
include("scripts/interactivity.jl")
include("scripts/plotting.jl")

# Display the resulting figure in its own window
display(figure);
