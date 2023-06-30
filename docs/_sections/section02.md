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
The *Raster File* is a file containing the dimensions of the study area, and is commonly one of the covariate data rasters. This raster will be a file on the device, and there are options to add the raster file from the device, remove the raster file, and export the template raster file. The models will be clipped to the extent of this raster. 

### Number of Multiprocessing Tiles (optional)
This is an optional input. Multiprocessing tiles provide a way to divide up a model's processing so that it can run faster and more efficiently. For more information about multiprocessing, see [SyncroSim Reference - Multiprocessing](https://docs.syncrosim.com/how_to_guides/modelrun_multiproc.html)
<br>

<p id="heading02"> <h2>Covariate Data</h2> </p>

The **Covariate Data** datasheet contains spatial information on each covariates. This information can be input directly, or imported from a csv (right click the datafeed -> Import) with the headers "CovariatesID"(covariate name), "RasterFilePath", "ResampleMethod" (optional), and "AggregationMethod" (optional).

### Covariate
In the *Covariate* column, the names of covariates that will be used in the *Scenario* are listed. By default, this column is empty. 
### Raster File
The *Raster File* column shows the file names of each covariate raster. 
### Resample Method
This *Resample Method* column defines which resampling method each raster will be using. This column can be created by right-clicking a covariate or raster file and selecting "Resample Method". Resampling options include Nearest Neighbor, Bilinear, Cubic, Cubic Spline, and Lanczos. To learn more about resampling rasters, see [PDF - Resampling Methods, MicroImages.com](https://www.microimages.com/documentation/TechGuides/77resampling.pdf)
### Aggregation Method
The *Aggregation Method* column defines which aggregation method each raster will be using. Similarly to *Resample Method*, this column can be created by right-clicking a covariate or raster file and selecting "Aggregation Method". Options for aggregation methods include Mean, Max, Min, and Majority. The aggregation method relates to changing the resolution of the raster. For more information about these aggregation methods see [Summarize Raster Within (Raster Analysis) - ESRI](https://pro.arcgis.com/en/pro-app/latest/tool-reference/raster-analysis/summarize-raster-within.htm) and [Aggregate (Spatial Analyst)](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/aggregate.htm).
<br>

<p id="heading03"> <h2>Field Data</h2> </p>

The **Field Data** datasheet contains information about species presence and absence locations.

### X
The *X* column represents X coordinates in the coordinate system defined in the **Options** tab. 
### Y
The *Y* column represents Y coordinates in the coordinate system defined in the **Options** tab. 
### Response
The *Response* column identifies whether there was an absence or presence detected at each site. 1 represents a presence detected, while 0 represents an absence detected. 
### Site
The *Site* column represents the site ID number. The site column can be added by right clicking a cell in the table. 
### Use In Model Evaluation
*Use In Model Evaluation* is an option that can be chosen for each row in this column, and can be added to the **Field Data** table by right clicking a row in the table and selecting "Use In Model Evaluation". Selecting "Yes" will save certain points for testing the *Scenario's* models, and selecting "No" will allow certain points to be used in the modeling process. 
### Model Selection Split
The *Model Selection Split* column reserves selected data from the model fitting process, but these data will still be included in model evaluation and the split data are stratified by the response. 
### Weight
The *Weight* column allows for weighting of certain points differently in the model. This column can be added to the table by right clicking a row in the table and selecting "Weight". 
<br>

<p id="heading04"> <h2>Options</h2> </p>

The **Options** datasheet controls some of the *Scenario's* settings relating to the **Field Data**.

### Authority Code (e.g., EPSG:4326)
The *Authority Code* contains coordinate system information for the scenario. This entry represents the coordinate system that the field data are in. 
### Aggregate or Weight Data
The *Aggregate or Weight Data* argument specifies whether the **Field Data** should be aggregated or weighted.
### Generate background sites
This argument specifies whether background sites should be generated within **SyncroSim** for the *Scenario*. 
### Number of background sites added to field data
If "Yes" is selected in the *Generate background sites* argument, the *Number of background sites added to field data* will specify the number of background points that should be generated for this *Scenario*. 
### Method used for background site generation
The *Method used for background site generation* will either be Kernel Density Estimate (KDE) or Minimum Convex Polygon (MCP). 
### KDE background surface method
The *KDE background surface method* will either be continuous or binary if the *Method used for background site generation* is KDE. 
### Isopleth threshold used for binary mask creation (%)
The *Isopleth threshold used for binary mask creation (%)* creates a binary mask using an isopleth threshold specified by the user if the KDE background surface method is binary or if the method used for background site generation is MCP. 

<p id="heading05"> <h2>Validation Options</h2> </p>

The **Validation Options** datasheet sets parameters for model validation.

### Split data for model training and testing
Checking "Yes" in the *Split data for model training and testing* tab will allow some of the data to be withheld from modeling and used for validating the models.
### Proportion of data used for model training
The *Proportion of data used for model training* option accepts proportions as decimals. For example, if this proportion is 0.3, 30% of the *Field Data* will be withheld for model validation. 
### Use cross validation for model selection
Checking "Yes" in the *Use cross validation for model selection* tab will further split data into cross validation folds, testing the models on a subset of the data that was not used to train the model. 
### Stratify cross-validation folds by the response
Checking "Yes" in the *Stratify cross-validation folds by the response* tab will ensure that the cross-validation folds will be representative of the distribution of the response values. 
### Number of cross-validation folds
The *Number of cross-validation folds* tab specifies how many cross-validation folds will be included in the *Scenario*.
<p id="heading06"> <h2>Spatial Multiprocessing</h2> </p>

The **Spatial Multiprocessing** datasheet contains a tiling raster used for spatial multiprocessing. 
### Tiling raster
The *Tiling raster* accepts a raster that represents the tiles that the data can be split into for parallel processing. 