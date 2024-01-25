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
	),
	Legend=(
		framevisible=false,
		gridshalign=:left,
		padding=(0.0f0, 0.0f0, 0.0f0, 0.0f0),
		labelcolor="#8f8f8f",
		titlecolor="#8f8f8f",
	),
	Colorbar=(tickcolor="#8f8f8f",),
	size=(960, 600),
)
