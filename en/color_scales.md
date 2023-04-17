# Color palettes
*Started on 28.3.2023 by Emmi. Feel free to contribute and add info and other resources here!*

Color palettes can be divided to sequential, diverging, cyclic, and qualitative. It's usually best to only use the type of palettes that match the type of the data, e.g. using diverging and cyclic palettes can look confusing if your color scale should represent percentages. Qualitative color palette's on the other hand usually provide better distinction between the colors than other types. Matplotlib (package for Python) has a nice introduction to these different kinds of colormaps in [Choosing colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html). 

It can also be good to consider if people with color vision deficiency (CVD) can interpret your figures. In northern Europe, 1 in 12 men and 1 in 200 women have a colour vision deficiency, see short intro into the topic [here](https://www.nature.com/articles/d41586-021-02696-z). Top two tips is not to use rainbows and to pair red and green (altough might be ok if they have different hues)

## Sequential

In general, sequential color scales were the lightness value increases monotonically through the colormaps and the change of color is perceptually uniform, are easier to comprehend than e.g. the infamous rainbow. With this kind of plots you're also ok if the figure should be readable also in greyscale. Sequential palettes can be also used as binned palettes, when just a finite number of colors (usually with equal spacing) are selected.

![sequential](https://matplotlib.org/stable/_images/sphx_glr_colormaps_015_2_0x.png)

### Viridis
The above figure is of colormaps in Matplotlib, but the [viridis](https://sjmgarnier.github.io/viridis/) package in R has the same colormaps and a few more. All colormaps in the viridis package should be interpretable for also people with the most common CVDs.

### Scico 
[scico](https://github.com/thomasp85/scico) has the same idea as Viridis in that all of its color palettes should be perceptually uniform, also for people with common CVDs. Scico has a bit more options for sequential palettes, and also diverging and cyclic palettes.

### Colorbrewer
[Colorbrewer](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3) has several sequential, diverging, and qualitative colormaps and many of these are also colorblind friendly. There's also a [colorbrewer R package](https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html).

## Diverging
As already mentioned, both Scico and Colorbrewer provide diverging color palettes. The R package for scico has a useful built-in option for setting the midpoint value without additional computations. 

## Qualitative

### Colorbrewer
As mentioned, colorbrewer has also qualitative colormaps.

### Matplotlib color palettes for 20+ colors
Using 20 or more colors may not be the neatest option for final figures, but can be sometimes useful. In Python Matplotlib has several good qualitative palettes ([Choosing colormaps in Matplotlib](https://matplotlib.org/stable/tutorials/colors/colormaps.html)). Especially if you need lots of colors, the tab20 is the best one I've found so far (especially if you happen to have paired data). tab20b and tab20c provide other options and can be partially combined for more colors. You can copy the palettes to R from below:

![qmatplotlib](https://matplotlib.org/stable/_images/sphx_glr_colormap_reference_006_2_0x.png)

```
tab20 <- c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78",
		   "#2ca02c","#98df8a","#d62728","#ff9896",
		   "#9467bd","#c5b0d5","#8c564b","#c49c94",
		   "#e377c2","#f7b6d2","#7f7f7f","#c7c7c7",
		   "#bcbd22","#dbdb8d","#17becf","#9edae5")

tab20b <- c("#393b79", "#5254a3", "#6b6ecf", "#9c9ede",
			"#637939", "#8ca252", "#b5cf6b", "#cedb9c",
			"#8c6d31", "#bd9e39", "#e7ba52", "#e7cb94", 
			"#843c39", "#ad494a", "#d6616b", "#e7969c", 
			"#7b4173", "#a55194", "#ce6dbd", "#de9ed6")

tab20c <- c("#3182bd", "#6baed6", "#9ecae1", "#c6dbef", 
			"#e6550d", "#fd8d3c", "#fdae6b", "#fdd0a2", 
			"#31a354", "#74c476", "#a1d99b", "#c7e9c0", 
			"#756bb1", "#9e9ac8", "#bcbddc", "#dadaeb", 
			"#636363", "#969696", "#bdbdbd", "#d9d9d9")
```

### Coolors
I've found that [Coolors](https://coolors.co/palettes/trending) has some nice looking color palettes you can browse. If you find something you like, just click on the three dots and choose export, to get e.g. an array of hex-codes you can use in R or Python.

### Okabe-Ito color palette
Okabe-Ito color palette was dewigned to be accessible for people with people with common CVDs. Check the [R package](https://easystats.github.io/see/reference/scale_color_okabeito.html) or copy the list of hex codes from below. 

![okabe-ito](https://jfly.uni-koeln.de/color/image/pallete.jpg)

```
okabeito <- c("#000000","#E69F00","#56B4E9","#009E73",
              "#F0E442","#0072B2","#D55E00","#CC79A7")
```

## Check readability for people with color vision deficiency

In R you can install e.g. [colorblindcheck](https://cran.r-project.org/web/packages/colorblindcheck/vignettes/intro-to-colorblindcheck.html) package to see how your color palette looks like to someone with Deuteranopia, protanopia, or tritanopia. If execute the code below, you'll see that `tab20` is not that good for people with CVDs (but color palettes with this many colors rarely are).

```
install.packages("colorblindcheck")
library(colorblindcheck)
palette_check(tab20), plot=TRUE)
```

If you'd like to check how your figure looks like, you can use
[Coblis](https://www.color-blindness.com/coblis-color-blindness-simulator/). Just copy-paste your figure there and it shows how it looks to people with CVDs.

## Summary of R packages (mentioned above) for color palettes
* [viridis](https://sjmgarnier.github.io/viridis/)
* [colorbrewer](https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html)
* [scico](https://github.com/thomasp85/scico)
* [scale\_color\_okabeito](https://easystats.github.io/see/reference/scale_color_okabeito.html)
* [colorblindcheck](https://cran.r-project.org/web/packages/colorblindcheck/vignettes/intro-to-colorblindcheck.html)

## Using color scales and other plotting information for R
* [ggplo2-book](https://ggplot2-book.org/scale-colour.html) for basics and general grammar for plotting with ggplot2.
	* Chapter 11 (Colour scales and legends) gives a little colour theory and describes e.g. how to apply different types of colour scales
* [R graphics Cookbook](https://r-graphics.org) for more than 150 receips for various types of graphs
