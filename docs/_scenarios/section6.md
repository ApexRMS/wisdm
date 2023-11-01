---
layout: default
section: 6
title: Output Options
permalink: reference/output-options
---

<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    <li>Library</li>
    {% for section in site.library %}
        {% if section.section != page.section %}
            <a href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% else %}
            <a class="selected" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% endif %}
    {% endfor %}
    <li>Project</li>
    {% for section in site.project %}
        {% if section.section != page.section %}
            <a href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% else %}
            <a class="selected" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% endif %}
    {% endfor %}
    <li>Scenario</li>
    {% for section in site.scenarios %}
        {% if section.section != page.section %}
            <a class="indent1" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% else %}
            <a class="selected" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% endif %}
    {% endfor %}
</div>

# **Output Options**

The **Output Options** tab is a *Scenario Datasheet*.

In the SyncroSim UI, the **Output Options** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

For more information about viewing outputs, see the [SyncroSim documentation](https://docs.syncrosim.com/how_to_guides/results_overview.html).

<br>

### **Make Probability Map**
Determines whether to create a probability map with the extent of the template raster for the **Scenario** ("Yes"). The map indicates the probability of finding a given species in each pixel/cell, given the constraints and biases of the sampling design. 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Make Residuals Map**
Determines whether to create a model deviance residuals map for the **Scenario** and inform the user about any issues with the model fit, should there be a spatial pattern in residuals ("Yes").

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make Multivariate Environmental Similarity Surface (MESS) Map**
Determines whether to create a MESS map for the **Scenario** that shows how well each point fits into the ranges of the points to which the model was fit ("Yes"). Negative values in this map indicate points outside of the training data ranges.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make Most Dissimilar Variable (MoD) Map**
Determines whether to create a MoD map for the **Scenario** and indicate which covariate was the furthest from the range of observations used for model training ("Yes").

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>