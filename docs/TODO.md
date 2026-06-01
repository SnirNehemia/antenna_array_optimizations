
# 27-05-2026

* level support from config

I want to make sure of several things:
1. when calculating the mean in _directive_cost in @src/cost/cost_function.py , is the mean claculated on the linear scale pattern? (since it's physically wrong to calculate mean straight on logarithmic units)
2. add a easy way in the config to switch to a cost function that optimize over the min/max value of the directive in the angular window (in addition to the current "mean" method)
3. where and how do you use level from 
4. does the "wrap around" of angles treated here? If I ask for a theta = 5 and theta_width = 20 (and phi=alpha), does it know to "wrap" around angle so that it optimize over theta=-5:0 at phi=alpha+180, and over theta=0:15 at phi=alpha? same question about phi? does it "wrap" at 0 or 360?

#

# 26/5/2026

* add to the manual_weights.py a dropdown for common tapering - for quick load comparison. - should I?

* get the complex array from the mail and demonstrate it wins. V

* add support for theta_width and phi_width V

* fix the tehta range in the raw results folder V

* in @manual_weights.py - in the metric section, add a table of max/min peak/null.... whatever ?

* move the metric to somewhere it is easily accessible. also, write a document about what we are doing here in the optimization and a separate one to show how we verify it.

* does the optimization "wrap around" in theta and phi? V

# 19-05-2026

* there is an error in the visualozation: the phi angle does not wrap around for the traget \phiW(\degs) angle. also, it should show a wrap arounf phi angles, which also is not shown now. V


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