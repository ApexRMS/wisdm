---
layout: default
section: 04
title: Models
permalink: section04
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;GLM</a>
            <a href="#heading02"> &emsp;Random Forest</a>
            <a href="#heading03"> &emsp;Maxent</a>
            <a href="#heading04"> &emsp;Model Outputs</a>
        {% endif %}
    {% endfor %}
</div>

# **Models**

This is the page for the **Models** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

<br>
<br>

<p id="heading01"> <h2>GLM</h2> </p>

The **GLM** datasheet contains information about the .

### Select Best Predictors

### Simplification Method

### Consider Squared Terms

### Consider Interactions

<br>

<p id="heading02"> <h2>Random Forest</h2> </p>

The **Random Forest** datasheet contains information on.

### Evaluate covariate importance

### Calculate casewise importance

### Number of variables sampled at split

### Maximum number of nodes

### Number of trees

### Node size

### Normalize votes

### Calculate proximity

### Sample with replacement


<br>
<p id="heading03"> <h2>Maxent</h2> </p>

The **Maxent** datasheet contains information on.

### Memory allocation (GB)

### Maximum number of background points

### Number of processing threads

### View maxent user interface during run

### Save maxent formated input/output files

<br>

<p id="heading04"> <h2>Model Outputs</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### Model RDS

<br>

