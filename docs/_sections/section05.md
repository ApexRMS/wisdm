---
layout: default
section: 05
title: Output Options
permalink: section05
---

<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
    {% endfor %}
</div>

# **Output Options**

This is the page for the **Output Options** tab with content organized into headings and, optionally, subsections. For more information about viewing outputs, see [SyncroSim - Analyzing Results: Overview](https://docs.syncrosim.com/how_to_guides/results_overview.html).

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

### Make Probability Map
Selecting "Yes" on the *Make Probability Map* argument will create a probability map for the **Scenario** indicating the probability of finding the species in an area, given the constraints and biases of the sampling design.
### Make Residuals Map
Selecting "Yes" on the *Make Residuals Map* argument will create a model deviance residuals map for the **Scenario** and can inform the user about any issues with the model fit if there is a spatial pattern in residuals.
### Make Multivariate Environmental Similarity Surface (MESS) Map
Selecting "Yes" on the *Make Multivariate Environmental Similarity Surface (MESS) Map* argument will create a MESS map for the **Scenario** that shows how well each point fits into the ranges of the points to which the model was fit. Negative values in this map indicate points outside of the training data ranges. 
### Make Most Dissimilar Variable (MoD) Map
Selecting "Yes" on the *Make Most Dissimilar Variable (MoD) Map* argument will create a MoD map for the **Scenario**, and will indicate which covariate was the furthest from the range of observations used for model training.