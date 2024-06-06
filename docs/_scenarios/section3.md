---
layout: default
section: 3
title: Data Preparation
permalink: reference/data-preparation
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
            <a href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% else %}
            <a class="selected" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
            <a href="#heading01"> &emsp;&emsp;&emsp;Template Raster</a>
            <a href="#heading02"> &emsp;&emsp;&emsp;Covariate Data</a>
            <a href="#heading03"> &emsp;&emsp;&emsp;Field Data</a>
            <a href="#heading04"> &emsp;&emsp;&emsp;Field Data Options</a>
            <a href="#heading05"> &emsp;&emsp;&emsp;Background Data Options</a>
            <a href="#heading06"> &emsp;&emsp;&emsp;Validation Options</a>
            <a href="#heading07"> &emsp;&emsp;&emsp;Spatial Multiprocessing</a>
        {% endif %}
    {% endfor %}
</div>

# **Data Preparation**

The **Data Preparation** tab groups the following *Scenario Datasheets*:
* Spatial Data: 
    * Template Raster
    * Covariate Data
    * Restriction Raster
* Field Data
    * Field Data Options
* Validation Options
* Spatial Multiprocessing

In the SyncroSim UI, the **Data Preparation** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

<br>

<p id="heading01"> <h2><b>Template Raster</b></h2> </p>

The **Template Raster** *Datasheet* can be found under the **Spatial Data** tab and contains information about the dimensions of the study area.

### **Raster File**
Defines the template raster (e.g., GEOTiff) file containing the target dimensions of the study area. The *Covariate Raster Files* will be clipped to the extent of this raster and will have the same spatial reference and resolution (cell size) as this raster. In the SyncroSim UI, the Raster File can be selected by clicking on the yellow folder icon, removed by clicking on the X icon, and exported by clicking on upwards arrow icon. 

### **Number of Multiprocessing Tiles**
*Optional*. Sets the number of multiprocessing tiles to divide up a model's processing across multiple cores. This option allows to model to run faster and more efficiently. For more information about multiprocessing, see [SyncroSim Reference - Multiprocessing](https://docs.syncrosim.com/how_to_guides/modelrun_multiproc.html).

<br>

<p id="heading02"> <h2><b>Covariate Data</b></h2> </p>

The **Covariate Data** *Datasheet* can be found under the **Spatial Data** tab and contains the spatial information about each covariate. Covariates, also known as predictors, represent spatial and environmental data (*e.g.*, precipitation and slope) that models rely on to infer spatial patterns to predict habitat suitability for a species. 

In the SyncroSim UI, this information can be input directly or imported from a CSV file (right-click the **Covariate Data** page and select Import) that contains the four columns described below. If importing the data into SyncroSim, it must contain the following headers and corresponding information: "CovariatesID" (see *Covariate*), "RasterFilePath" (see *Raster File*), "ResampleMethod" (see *Resample Method*), and "AggregationMethod" (see *Aggregation Method*). 

### **Covariate**
Defines the names of the covariates that will be used in the *Scenario*. 

### **Raster File**
Defines the path to the covariate's raster (.tif) file. In the SyncroSim UI, the Raster File can be selected by clicking on the yellow folder icon, removed by clicking on the X icon, and exported by clicking on upwards arrow icon. 

### **Resample Method**
*Optional*. Defines which resampling method each raster will use. Resampling interpolates cell values when transforming the covariate raster to the coordinate reference system (CRS) or resolution of the template raster. Resampling options include: Nearest Neighbor, Bilinear, Cubic, Cubic Spline, and Lanczos. To learn more about resampling rasters, see [Resampling Methods](https://www.microimages.com/documentation/TechGuides/77resampling.pdf). In the SyncroSim UI, this column is initially hidden, and can be made visible by right-clicking anywhere inside the **Covariate Data** *Datasheet* and selecting "Resample Method" from the context menu. 

<div class=indentation> 
    <i>Default:</i> If a covariate is categorical, "Nearest Neighbor" is the default resampling method. If a covariate is not categorical, "Bilinear" is the default resampling method.
</div>

### **Aggregation Method**
*Optional*. Defines which aggregation method each raster will be using. Aggregation of covariate layers occurs when the resolution of the raster needs to be decreased to match the template raster. Aggregating ensures that the data within each cell is preserved when decreasing the resolution of the covariate raster to match the template raster (for example, re-scaling a 10 m resolution covariate raster to a 100 m resolution template raster). Options for aggregation methods include Mean, Max, Min, and Majority. For more information about these aggregation methods, see [Summarizing Raster](https://pro.arcgis.com/en/pro-app/latest/tool-reference/raster-analysis/summarize-raster-within.htm) and [Aggregating Rasters](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/aggregate.htm). In the SyncroSim UI, this column is initially hidden, and can be made visible by right-clicking anywhere inside the **Covariate Data** *Datasheet* and selecting "Aggregation Method" from the context menu.

<div class=indentation> 
    <i>Default:</i> If a covariate is categorical, "Majority" is the default aggregation method. If a covariate is not categorical, "Mean" is the default aggregation method. 
</div>

<br>

<p id="heading03"> <h2><b>Restriction Raster</b></h2> </p>

The **Restriction Raster** *Datasheet* can be found under the **Spatial Data** tab and contains information pertaining to the optional restriction layer.

### **Raster File**
This is an optional input that accepts a raster (.tif) file that gets multiplied across the probability raster during the apply model stage. It is recommended that this is a binary raster, where 0 indicates restricted areas where prediction is not necessary. The *Raster File* could otherwise be a probability raster (i.e., values ranging from 0.0 to 1.0), which would act as a weighted restriction layer, adjusting the values of all pixels based on the pixel restriction values.

### **Description**
This is also an optional input that can be used to detail the corresponding **Raster File**.

<br>

<p id="heading04"> <h2><b>Field Data</b></h2> </p>

The **Field Data** *Datasheet* can be found under the **Field Data** tab and contains information about species presence and absence locations.

### **X**
Defines the X coordinates (longitude) in the coordinate system defined in the <a href="#heading05">Options</a> tab. 

### **Y**
Defines the Y coordinates (latitude) in the coordinate system defined in the <a href="#heading05">Options</a> tab. 

### **Response**
Defines whether there was an absence or presence detected at a given site, where 1 represents a location where the species was present and 0 represents a location where the species was absent. 

### **Site**
*Optional*. Sets the site ID number. In the SyncroSim UI, this column is initially hidden, and can be made visible by right-clicking anywhere inside the **Field Data** *Datasheet* and selecting "Site" from the context menu. 

### **Use In Model Evaluation**
*Optional*. Determines whether a given row of data was saved for testing the Scenario's models (Yes) or used in the model training process (No). This column is auto-filled based on the parameters defined under <a href="#heading07">Validation Options</a>.

### **Model Selection Split**
*Optional*. Sets in which cross-validation fold a given row of data will be included. This column is auto-filled based on the parameters defined under <a href="#heading07">Validation Options</a>. Its values will range from 1 to the defined *Number of cross-validation folds*. Rows of data marked for *Model Selection Split* will be reserved from the model training process but will still be included in model evaluation. Additionally, the split data can be stratified by the response variable to ensure the proportion of presence and absence points is balanced between the testing and training split provided *Stratify cross-validation folds by the response* is set to "Yes".

### **Weight**
*Optional*. Sets the weight of a given data point. This column is only used by the model if *Aggregate or Weight Data* in the <a href="#heading05">Options</a> *Datasheet* is set to "Weights". If multiple species presence/absence points are included in a single pixel, these points can be down-weighted so that pixels with multiple points can have an equal weight to pixels with a single point. In the SyncroSim UI, this column is initially hidden, and can be made visible by right-clicking anywhere inside the **Field Data** *Datasheet* and selecting "Weight" from the context menu. 

<br>

<p id="heading05"> <h2><b>Field Data Options</b></h2> </p>

The **Field Data Options** *Datasheet* can be found under the **Field Data** tab and controls some of the *Scenario's* settings relating to the **Field Data**.

### **Authority Code (e.g., EPSG:4326)**
*Optional*. Sets the coordinate reference system (CRS) information for the **Field Data** points. This entry represents the CRS that the field data are in. If left blank, **WISDM** assumes that the field data are in the same CRS as that of the template raster. 

### **Aggregate or Weight Data**
*Optional*. Determines whether the **Field Data** points should be "Aggregated" or "Weighted" in the event that multiple points fall within the same pixel. Aggregating will combine multiple points into one, while weighting will ensure that all of the points are downweighted so that combined, they have the same influence as just one point. 

<br>

<p id="heading06"> <h2><b>Background Data Options</b></h2> </p>

The **Background Data Options** *Datasheet* controls some of the *Scenario's* settings relating to the **Field Data**.

### **Generate background sites**
Defines whether background sites should be generated for the *Scenario*. Background sites are also referred to as pseudo-absences and represent absences of a species. This information allows the models to compare environmental spaces where a species can and cannot be found. Background sites are often used when true absence data are not available for a species.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Number of background sites added to field data**
*Optional*. Sets the number of background points that should be generated for this *Scenario*. Depends on *Generate background sites* being set to "Yes". 

<div class=indentation> 
    <i>Default:</i> Sum of presence points from the <b>Field Data</b> <i>Datasheet</i>.
</div>

### **Method used for background site generation**
*Optional*. Defines the background surface method for generating background points based on the species presence values. Two options are available: "Kernel Density Estimate (KDE)" or "Minimum Convex Polygon (MCP)".

<div class=indentation> 
    <i>Default:</i> Kernel Density Estimate (KDE).
</div>

### **KDE background surface method**
*Optional*. Defines the type of *KDE background surface method* to be used if the *Method used for background site generation* is set to "Kernel Density Estimate (KDE)". Two options are available: "Continuous" or "Binary" . 

<div class=indentation> 
    <i>Default:</i> Continuous.
</div>

### **Isopleth threshold used for binary mask creation (%)**
*Optional*. Sets the threshold value that is required to be used when creating a binary mask if the *Method used for background site generation* is set to "Kernel Density Estimate (KDE)" and *KDE background surface method* is set to "Binary", or if *Method used for background site generation* is set to "Minimum Convex Polygon (MCP)". The isopleth threshold value is expressed as an integer, so an isopleth of 99% would require an input of 99 for this argument. 

<div class=indentation> 
    <i>Default:</i> 95.
</div>

<br>

<p id="heading07"> <h2><b>Validation Options</b></h2> </p>

The **Validation Options** *Datasheet* sets parameters for model validation. In this documentation, model processing, model fitting, and model training are used interchangeably and refer to modeling using training data. Model validation and model testing are used interchangeably and refer to the evaluation of model predictive performance using testing data.

### **Split data for model training and testing**
Determines whether to withhold data from modeling to be used for model validation.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Proportion of data used for model training**
Sets the proportion of the *Field Data* to be used for model training if *Split data for model training and testing* is set to "Yes". Proportions should be entered decimals. For example, setting this parameter to 0.7, 70% of the *Field Data* will be used for model training and 30% of the *Field Data* will be withheld for model validation (testing). 

<div class=indentation> 
    <i>Default:</i> 0.5.
</div>

### **Use cross validation for model selection**
Determines whether to split data into cross validation folds so that models are tested on a subset of the data that was not used for model training.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Stratify cross-validation folds by the response**
Determines whether the cross-validation folds should be representative of the distribution of the response values. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Number of cross-validation folds**
Sets the number of cross-validation folds to include in the *Scenario*. 

<div class=indentation> 
    <i>Default:</i> 10.
</div>

<br>

<p id="heading08"> <h2><b>Spatial Multiprocessing</b></h2> </p>

The **Spatial Multiprocessing** *Datasheet* contains a tiling raster used for spatial multiprocessing.

### **Tiling raster**
Defines the raster (e.g., GEOTiff) file representing the tiles the data can be split into for parallel processing. 

<br>