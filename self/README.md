# Band Contrast Sim

This program simulates band contrast images based on AFM data. Once simulated, a mapping between the measured band contrast and the simulated band contrast is created.

## To Compile

Run the makefile with `make` to compile the program.

## Before Running

This program depends on a specific file structure. The following directories must exist:

`~/Facets/samples/{sample name}/Input`  
`~/Facets/samples/{sample name}/Output`

This program will attempt to read AFM, band contrast, and (possibly) facet data from these directories.

AFM:  
`~/Facets/samples/{sample name}/Input/{sample name}.txt`  
Band Contrast:  
`~/Facets/samples/{sample name}/Output/BandContrast.txt`  
Facets:  
`~/Facets/samples/{sample name}/Output/thetamap_0.001.txt`  
`~/Facets/samples/{sample name}/Output/phimap_0.001.txt`

## To Run

The executable takes 6 arguments: SampleName, tiltAngleDeg, lightAngleDeg, widthInUM, alpha (fraction of light for linCombo [0-1]), readFacets (0 to calculate facets, 1 to read from files)

An example run command is as follows:  
`./BandContrastSim EC2_000 70 20 8 0.8 1`

That is:  
A 70 degree sample tilt;  
the light shines at 20 degrees (70 + 20 = 90, so the light is shining straight down in the tilted view);  
the sample is 8um wide;  
the linear combination of light and detector views is 80% light and 20% detector;  
the facets will be read from previously made facet files; 

If the facet files do not exist, enter a 0 for the last argument and they will be created (longer startup time). This only need to be done once. Facet files can then be read in to save time on subsequent trials.

## Mapping

### Manual Mapping

The mapping involves 12 parameters: a0-a5 and b0-b5. Let x' and y' be AFM coordinates, x and y be measured band contrast coordinates, and X and Y be the average x and y respectively. The 12 parameters are used to create a mapping as follows:

`x' = a0 + a1(x - X) + a2(y - Y) + a3(x - X)^2 + a4(y - Y)^2 + a5(x - X)(y - Y)`  
`y' = b0 + b1(x - X) + b2(y - Y) + b3(x - X)^2 + b4(y - Y)^2 + b5(x - X)(y - Y)`

When prompted, you will be able to enter values for these parameters. Enter space separated numbers to update the parameters in order from a0-a5-b0-b5.  
If you do not want to change a value from the previous run, enter a '.'.  
You need only enter as many values up until you no longer want to change something. That is, If you only want to change a1, you may enter:

`. {a1's new value}`

If you want to change, say a1 and b0, you may enter something like:

`. {a1's new value} . . . . {b0's new value}` 

### Amoeba

Once you have an output image that is visibly close to the desired mapping, you can run amoeba. This will attempt to find the best local by minimizing chi squared. To run amoeba, use the `a {fitLevel}` command. The `fitLevel` corresponds to the number of parameters you wish amoeba to try and fit at once. Possible amoeba commands are as follows:

`a 1` -> Fit only the position parameters  
`a 2` -> Fit all previous + scale parameters  
`a 3` -> Fit all previous + skew parameters  
`a 4` -> Fit all previous + stretch parameters  
`a 5` -> Fit all previous + curve parameters  
`a 6` -> Fit all previous + splay(?) parameters

## The Output

You can find the output in: `~/Facets/samples/{sample name}/simOutput/`.

You will find several greyscale images (.pgm). When mapping, openning `overlap.pgm` will show the most recent mapping of the measured band contrast onto the simulated band contrast. `coeffNChi.txt` contains the coefficients you have tried to use for mapping.  
**NOTE**: the last number on a line is the chi squared for those coefficients. Ensure you do not copy and past it into the program while mapping.
