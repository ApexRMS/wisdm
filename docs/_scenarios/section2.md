---
layout: default
section: 2
title: General
permalink: reference/general
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
            <a href="#heading01"> &emsp;&emsp;&emsp;Summary</a>
            <a href="#heading02"> &emsp;&emsp;&emsp;Pipeline</a>
            <a href="#heading03"> &emsp;&emsp;&emsp;Datafeeds</a>
        {% endif %}
    {% endfor %}
</div>

# **General**

The **General** tab groups the following *Scenario Datasheets*::
* Summary
* Pipeline
* Datafeeds

In the SyncroSim UI, the **General** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

<br>

<p id="heading01"> <h2><b>Summary</b></h2> </p>

The **Summary** *Datasheet* contains general information about the *Scenario*. 

### **Name**
Defines the name of the *Scenario*. The name of the *Scenario* can be defined when creating a new *Scenario* and modified by editing this field. In the SyncroSim UI, the "Name" can also be edited by right-clicking the *Scenario* on the *Library Explorer* window and selecting ***Rename...***.

### **Owner**
Defines the owner of the *Scenario*.

### **Description**
Provides a description of the *Scenario*, including objectives, species names, geographic extent, and other pertinent information to be captured.

### **Library**
Defines the *Library* to which the *Scenario* belongs. The *Library* is the highest level of organization within the *Library Explorer* window, and contains the *Project*. In the SyncroSim UI, the *Library* name can be modified by right-clicking the *Library* on the *Library Explorer* window, selecting ***Properties*** from the context menu, and editing the "Name" field .
> R and Python configurations are controlled at the *Library* level.

### **Project**
Defines the *Project* to which the *Scenario* belongs. The *Project* is the second highest level of organization within the *Library Explorer* window. It is nested within the *Library* level and contains all *Scenarios* for a specific project. In the SyncroSim UI, the default *Project* name "Definitions" can be modified by right-clicking the *Project* on the *Library Explorer* window, selecting ***Properties*** from the context menu, and editing the "Name" field.
> Covariates are defined at the *Project* scope.

### **Last modified**
Defines when the *Scenario* was last modified. It is presented in the following format: Month-Day-Year Hour:Minute:Second AM/PM. 

### **Auto generation tags**
Used to autogenerate *Scenarios* for factorial combinations of parameter inputs.

### **Read only**
Checking "Read only" prevents the *Scenario* from being edited. 

### **Merge Dependencies**
Dependencies allow different steps of the pipeline to be run in different *Scenarios*. For more information about dependencies, see the [SyncroSim documentation](https://docs.syncrosim.com/how_to_guides/properties_dependencies.html). Checking "Merge Dependencies" allows the dependencies for the source *Scenario(s)* to be merged rather than being prioritized by the order they appear in the list of dependencies.

### **Ignore Dependencies...**
Selecting "Ignore Dependencies" will lead to a panel showing *Datafeeds* and *Packages*. Checking the box "Ignore" allows the *Scenario* to ignore the selected *Package(s)* when running. 

<br>

<p id="heading02"> <h2><b>Pipeline</b></h2> </p>

The **Pipeline** *Datasheet* is a SyncroSim Core *Datasheet* that controls which model transformers to run and their run order. For more details, see [Scenario Pipeline](https://docs.syncrosim.com/reference/ds_scenario_pipeline.html). 

### **Stage**
Defines the transformers that will be run in the *Scenario*. The general stages in **WISDM** include:
1. Prepare Multiprocessing
2. Spatial Data Preparation
3. Site Data Preparation
4. Background Data Generation
5. Prepare Training/Testing Data
6. Variable Reduction
7. Models 
    * Boosted Regression Trees (BRTs)
    * Generalized Linear Model (GLM)
    * Maxent
    * Random Forest
8. Apply Model
9. Ensemble Model

Here are some short descriptions of each stage:
1. <u>Prepare Multiprocessing</u>: 
    * Generates a tiling raster that will be used by subsequent transformers to divide spatial data into smaller chunks for parallel processing. 
        * Required inputs: Template Raster File. 
        * Outputs: Spatial Multiprocessing Tiling Raster.

2. <u>Spatial Data Preparation</u>: 
    * Prepares Covariate and Restriction rasters to align with the extent, resolution, and reference system of the Template raster. 
        * Required inputs: Template Raster File, Covariate data. 
        * Outputs: Updated covariate data.

3. <u>Site Data Preparation</u>: 
    * Field Data are projected and clipped to the resolution and extent of the Template Raster. If Aggregate or Weight is selected under **Field Data Options**, the field data records are aggregated or weighted based on their spatial distribution. Covariate values are then extracted for the remaining field data records. 
        * Required inputs: Template Raster File, Field Data, Covariate Data. 
        * Outputs: Updated Field Data and Site Data.

4. <u>Background Data Generation</u>: 
    * Generates background sites (pseudoabsences), extracts covariate data for the pseudoabsence sites and updates the Field Data and Site Data to include pseudoabsence data. 
        * Required inputs: Template Raster File, Field Data, Covariate Data, Site Data. 
        * Outputs: Updated Field Data and updated Site Data. 

5. <u>Prepare Training/Testing Data</u>: 
    * Divides Field Data into training and testing sets and/or cross validation folds based on the provided Validation Options arguments. 
        * Required inputs: Field Data, Site Data. 
        * Outputs: Updated Field Data. 

6. <u>Variable Reduction</u>: 
    * Reduces correlated covariates through either an automated approach (using Variance Inflation Factors) or an interactive approach which involves user selection of covariates through a correlation viewer application. 
        * Required inputs: Field Data, Site Data. 
        * Outputs: Updated Retained Covariate List and a Covariate Selection Matrix if the interactive approach was used.

7. <u>Fit Models</u> (Boosted Regression Trees, Generalized Linear Model, Maxent, Random Forest): 
    * Fits the selected model and generates validation and diagnostic summaries and figures. 
        * Required inputs: Field Data, Validation Options, Site Data, Retained Covariate List. 
        * Outputs: Updated Model Outputs.

8. <u>Apply Model</u>: 
    * Generates map outputs for the selected model(s). Maps are generated based on the model, selected Output Options, and the provided Covariate Data. 
        * Required inputs: Template Raster File, Model Outputs, Covariate Data. 
        * Outputs: Map Outputs.

9. <u>Ensemble Model</u>: 
    * Generates ensembled outputs of the models run in the **Scenario**. Ensembles are generated based on the fitted models, apply models, and selected Ensemble Options. 
        * Required inputs: Model Outputs, Map Outputs. 
        * Outputs: Ensemble Outputs.

### **Run Order**
Sets the order in which the stages will be run.

### **Jobs**
Not currently used by **WISDM**. Within SyncroSim, it sets the maximum number of jobs per stage. Within **WISDM**, the number of jobs is defined by the number of multiprocessing tiles. 

<br>

<p id="heading03"> <h2><b>Datafeeds</b></h2> </p>

The **Datafeeds** *Datasheet* summarizes all the *Datafeeds* parameterized for the *Scenario*. 

### **Data**
Determines whether a given *Datasheet* has been parameterized. When a new *Scenario* is created, the "Data" field will be blank, as no data or information has been provided yet. Once a *Datasheet* <u>has</u> been parameterized, the "Data" field will be set to True and in the SyncroSim UI, a green checkmark will appear.

### **Package Name**
Defines to which *Package* a given *Datasheet* belongs. For **WISDM** models, three *Package* names will be displayed: "wisdm", "corestime", and "core".

### **View**
Displays the *Datasheet* and *Datafeed* names. In the SyncroSim UI, the name will be displayed as a hyperlink that leads to the respective *Datasheet*. 

### **Datafeed**
Defines the names of each *Datafeed* in the *Scenario*. In **WISDM**, the *Datasheets* and their respective *Datafeeds* include:
* General - Pipeline
    * Pipeline
* Data Preparation - Spatial Data
    * Template Raster
    * Covariate Data
    * Restriction Raster
* Data Preparation - Field Data
    * Field Data 
    * Field Data Options
    * Background Data Options
* Data Preparation - Validation Options
    * Validation Options
* Data Preparation - Spatial Multiprocessing
    * Spatial Multiprocessing
* Variable Reduction
    * Site Data
    * Retained Covariate List
    * Covariate Selection Options
* Models
    * GLM
    * Random Forest
    * Maxent
    * Boosted Regression Tree
    * Model Outputs
* Output Options
    * Output Options
    * Ensemble Options

### **Source Scenario**
Defines the *Scenario* from which each *Datasheet* is being drawn. 

<br>