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

### **Project**
Defines the *Project* to which the *Scenario* belongs. The *Project* is the second highest level of organization within the *Library Explorer* window. It is nested within the *Library* level, and contains all of the *Scenarios*. In the SyncroSim UI, the default *Project* name "Definitions" can be modified by right-clicking the *Project* on the *Library Explorer* window, selecting ***Properties*** from the context menu, and editing the "Name" field. 
> Covariates are defined at the *Project* scope.

### **Library**
Defines the *Library* to which the *Scenario* belongs. The *Library* is the highest level of organization within the *Library Explorer* window, and contains the *Project*. In the SyncroSim UI, the *Library* name can be modified by right-clicking the *Library* on the *Library Explorer* window, selecting ***Properties*** from the context menu, and editing the "Name" field .
> R and Python configurations are controlled at the *Library* level.

### **Last modified**
Defines when the *Scenario* was last modified. It is presented in the following format: Month-Day-Year Hour:Minute:Second AM/PM. 

### **Auto generation tags**
Used to autogenerate *Scenarios* for factorial combinations of parameter inputs.

### **Read only**
Checking "Read only" prevents the *Scenario* from being edited. 

### **Merge Dependencies**
Checking "Merge dependencies" allows the dependencies for the source *Scenario(s)* to be merged, rather than prioritized. Dependencies allow different steps of the pipeline to be run in different *Scenarios*. For more information about dependencies, see the [SyncroSim documentation](https://docs.syncrosim.com/how_to_guides/properties_dependencies.html).

### **Ignore Dependencies...**
Selecting "Ignore Dependencies" will lead to a panel showing *Datafeeds* and *Packages*. Checking the box "Ignore" allows the *Scenario* to ignore the selected *Package(s)* when running. 

<br>

<p id="heading02"> <h2><b>Pipeline</b></h2> </p>

The **Pipeline** *Datasheet* is a SyncroSim Core *Datasheet* that controls which model transformers to run and their run order. For more details, see [Scenario Pipeline](https://docs.syncrosim.com/reference/ds_scenario_pipeline.html). 

### **Stage**
Defines the transformers that will be run in the *Scenario*. The general stages in **WISDM** include:
1. Prepare Multiprocessing
2. Spatial Data Preparation
3. Data Preparation (Non-Spatial)
4. Variable Reduction
5. Models (Maxent, Random Forest, GLM)
6. Apply Model

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
Defines the names of each *Datafeed* in the *Scenario*. In **WISDM**, the *Datafeeds* and their respective *Datasheets* include:
* Data Preparation
    * Template Raster
    * Covariate Data
    * Field Data
    * Field Data - Options
    * Validation Options
* Variable Reduction
    * Site Data
    * Covariate Selection Options
    * Reduced Covariate List
* Models
    * GLM
    * Random Forest
    * Maxent
    * Model Outputs
* Output Options

### **Source Scenario**
Defines the *Scenario* from which each *Datasheet* is being drawn. 

<br>