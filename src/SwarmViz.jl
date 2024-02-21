module SwarmViz

using DataFrames
using Colors
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

	# run the scripts to set up the style and layout,
	# add interactivity, and create reactive plots
	#! Order of execution matters and should never be changed!
	include("scripts/plotstyle.jl")
	include("scripts/layout.jl")
	include("scripts/interactivity.jl")
	include("scripts/plotting.jl")

	# Display the resulting figure in its own window
	display(figure);
	return 0
end

end
