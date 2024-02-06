#This script sets up all elements of the layout

# set the title of the window that opens upon execution of the main script
GLMakie.activate!(; title="SwarmViz")

# create the figure and set up custom fonts for the UI
figure = Figure(;
	size=(960, 600),
	fonts=(;
		regular=joinpath(FONTFOLDER, "fira_sans_font", "FiraSans-Medium.ttf"),
		ui_font=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Regular.ttf"),
		ui_font_bold=joinpath(FONTFOLDER, "ubuntu_font", "Ubuntu-Bold.ttf"),
	),
)

# create the four high level areas of the layout
swarm_animation = Axis(
	figure[1, 1]; xlabel="X", ylabel="Z", autolimitaspect=1, alignmode=Outside()
)
metrics_grid = GridLayout(figure[1, 2]; alignmode=Outside())
time_grid = GridLayout(figure[2, 1:2])
controls = GridLayout(figure[3, 1:2])

#adjust the width of the columns
colsize!(figure.layout, 1, Relative(0.5))
colsize!(figure.layout, 2, Relative(0.5))

# TIMETEP

# create the sliders that controls the current timestep together
# with custom labels that allow changing the font
timestep = Slider(time_grid[1, 2]; range=timesteps, startvalue=1)
Label(time_grid[1, 1]; text="Timestep", halign=:left, font=:ui_font)
Label(time_grid[1, 3], @lift string($(timestep.value)); halign=:right, font=:ui_font)
#TODO add textbox to enter specific value?

# CONTROLS

# play/pause button with a reactive label
play_button_text = @lift $isplaying[] ? "PAUSE" : "PLAY"
play_button = Button(
	controls[1, 1]; label=play_button_text, width=60, height=60, font=:ui_font_bold
)

# sliders for fps and numnber of skipped frames
# with custom labels to allow changing the font
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
# make this column resize with the window but take up exactly 25% of the width
colsize!(controls, 2, Relative(0.25))

# toggles to display the convex hull, center and diameter of the swarm
animation_toggles = [Toggle(figure) for _ in 1:3]
controls[1, 3] = grid!(
	hcat(
		[
			Label(figure, l; halign=:left, font=:ui_font) for
			l in ["Convex Hull", "Center", "Diameter"]
		],
		animation_toggles,
	);
	default_rowgap=3,
	default_colgap=6,
	tellheight=true,
)

# toggles to control what is displayed on the robot markers in the animation
# as well as which threshold is applied to the hierarchical clustering for the coloring
robot_controls = GridLayout(controls[1, 4]; default_rowgap=3)
robot_controls_upper = GridLayout(robot_controls[1, 1]; default_rowgap=3)
robot_toggles = [Toggle(figure) for _ in 1:2]
robot_controls_upper[1, 1] = grid!(
	hcat(
		[
			Label(figure, l; halign=:left, font=:ui_font) for
			l in ["Collisions", "Clustering"]
		],
		robot_toggles,
	);
	default_rowgap=3,
	default_colgap=3,
	halign=:left,
)
manual_log_threshold = Textbox(
	robot_controls_upper[1, 2];
	placeholder="Enter Log Thresh.",
	# font=:ui_font,
	tellwidth=false,
	halign=:right,
    reset_on_defocus = true,
    validator = Float64,
)

robot_controls_lower = GridLayout(robot_controls[2, 1])
# the possible threshold values are adapted to the minimum and maximum in the current data
heightrange = @lift round.(
	range(log.(extrema(reduce(vcat, [c.heights for c in $data.clustering])))..., 3000),
	sigdigits=4,
)
log_threshold = Slider(
	robot_controls_lower[1, 2]; range=heightrange, startvalue=(@lift median($heightrange))
) #TODO: fixed choice of range after deciding on dissimilarity measure? Rethink range in general
Label(robot_controls_lower[1, 1], "Log Threshold"; halign=:left, font=:ui_font)
Label(
	robot_controls_lower[1, 3],
	@lift string($(log_threshold.value));
	halign=:right,
	font=:ui_font,
    tellwidth=true,
    width=30, #TODO fine tune depending on final range
)
# make this column resize with the window and take up all remaining available space
colsize!(controls, 4, Auto())

# buttons for IO
data_controls = GridLayout(controls[1, 5]; default_rowgap=6, default_colgap=6)
import_button, wall_button, collision_button, export_metrics_button =
	data_controls[1:2, 1:2] = [
		Button(figure; label=l, halign=:left, width=w, height=27, font=:ui_font) for
		(l, w) in zip(
			["Import Tracking", "Import Wall", "Import Collisions", "Export All"],
			[120, 90, 120, 90],
		)
	]

# METRIC PLOTS

# axes to hold the metric plots with a narrow one below them for the collision events
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
