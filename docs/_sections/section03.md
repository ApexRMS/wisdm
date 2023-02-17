---
layout: default
section: 03
title: Variable Reduction
permalink: section03
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;Site Data</a>
            <a href="#heading02"> &emsp;Covariate Selection Options</a>
            <a href="#heading03"> &emsp;Reduced Covariate List</a>
        {% endif %}
    {% endfor %}
</div>

# **Variable Reduction**

This is the page for the **Variable Reduction** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

<br>
<br>

<p id="heading01"> <h2>Site Data</h2> </p>

The **Template Raster** datasheet contains information about the dimensions of the study area.

### Site

### Covariate

### Value

<br>

<p id="heading02"> <h2>Covariate Selection Options</h2> </p>

The **Covariate Data** datasheet contains spatial information on each covariates.

### Display Highest Correlations

### Correlation Threshold

### Number of Plots


<br>

<p id="heading03"> <h2>Reduced Covariate List</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### Covariate

<br>

