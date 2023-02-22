---
layout: default
section: 02
title: Data Preparation
permalink: section02
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;Template Raster</a>
            <a href="#heading02"> &emsp;Covariate Data</a>
            <a href="#heading03"> &emsp;Field Data</a>
            <a href="#heading04"> &emsp;Options</a>
            <a href="#heading05"> &emsp;Validation Options</a>
            <a href="#heading06"> &emsp;Spatial Multiprocessing</a>
        {% endif %}
    {% endfor %}
</div>

# **Data Preparation**

This is the page for the **Data Preparation** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

<br>
<br>

<p id="heading01"> <h2>Template Raster</h2> </p>

The **Template Raster** datasheet contains information about the dimensions of the study area.

### Raster File

### Number of Multiprocessing Tiles (optional)

<br>

<p id="heading02"> <h2>Covariate Data</h2> </p>

The **Covariate Data** datasheet contains spatial information on each covariates.

### Covariate

### Raster File

### Resample Method

### Aggregation Method

<br>

<p id="heading03"> <h2>Field Data</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### X

### Y

### Response

### Site

### Use In Model Evaluation

### Model Selection Split

### Weight

<br>

<p id="heading04"> <h2>Options</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### Authority Code (e.g., EPSG:4326)

### Aggregate or Weight Data

<p id="heading05"> <h2>Validation Options</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### Split data for model training and testing

### Proportion of data used for model training

### Use cross validation for model selection

### Stratify cross-validation folds by the response

### Number of cross-validation folds

<p id="heading06"> <h2>Spatial Multiprocessing</h2> </p>

Some random text here just to show the structure of the page and how internal links work.

### Tiling raster