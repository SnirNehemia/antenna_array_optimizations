
# 19-05-2026

* there is an error in the visualozation: the phi angle does not wrap around for the traget \phiW(\degs) angle. also, it should show a wrap arounf theta angles, which also is not shown now.

* make it so when the cursor is on the plot, it displays the value there - the actual value is the one in diaply, depending on the polarization and mode (relative\absolute) the user chose.

# 15-05-2026

* Fix the Mercator projection we use here: maybe add an "importance" matrix that fixes the projection to account for that there are more pixels in the flat picture of the radiation pattern in the poles (each pixel there accounts for smaller ).
also, fix the rectangular preview of the goal on the radiation pattern to ccount for it (now it wouldn't be a pefect rectangle)

* Add a theta and phi angular width for the target for both scripts.

* Show the radiation pattern improvement with gradient for each initial value as gif - showing the radiation pattern, its target and the performance in an animated graph.



# 13-05-2026

add the option of crosspol

note to myself: it should be quite straight-forward.

* run new cst models: one of "clear" array, one obstracted - 

* check if co\crosspol is the one to go with - in Beast

* and one "no-array"