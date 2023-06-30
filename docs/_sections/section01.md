---
layout: default
section: 01
title: General
permalink: section01
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;Summary</a>
            <a href="#heading02"> &emsp;Pipeline</a>
            <a href="#heading03"> &emsp;Datafeeds</a>
        {% endif %}
    {% endfor %}
</div>
<br>

# **General**

This is the page for the **General** tab with content organized into headings and, optionally, subsections. The **General** tab can be found within the *Properties* pane by right-clicking the *Scenario*.

The **Sidebar Navigation Menu** lists all the headings within the page and all the Sections within the Reference guide. 
<br>
<br>

<p id="heading01"> <h2>Summary</h2> </p>

The **Summary** datasheet contains general information about the *Scenario*. 

### Name
This column details the name of the *Scenario*. The name of the *Scenario* can be defined when opening a new *Scenario* and typing in the File Name. It can be changed in this column, or by right-clicking the *Scenario* and choosing *Rename...*.
### Owner
The owner of the *Scenario* can be defined in this column. 
### Description
Here, the description of the *Scenario* can be summarized. Objectives, species names, geographic extent and other information can be included in this section. 

### Project
The *Project* is the second highest level of organization within the *Library Explorer*, beneath the *Library* level, and holds *Scenarios* under this level. This entry shows to which *Project* the current *Scenario* belongs. The default *Project* name is "Definitions", although this can be changed by right-clicking the *Project* and selecting "Properties". Additionally, covariates can be defined at the *Project* level. 

### Library
The *Library* is the highest level of organization within the *Library Explorer* and can hold *Projects*. R and Python configuration for each *Project* is controlled at the *Library* level. This entry shows to which *Library* the *Scenario* belongs. The *Project* name can be changed by right clicking the *Library* and selecting *Properties*. To learn more about *Libraries* and *Projects*, see [SyncroSim - Libraries, Projects & Scenarios: Overview](https://docs.syncrosim.com/how_to_guides/library_overview.html)
### Last modified
This entry shows when the *Scenario* was last modified, Month/Day/Year Hour:Minute:Second AM/PM. 
### Auto generation tags
*Auto generation tags* are not a requirement, but can help autogenerate *Scenarios* for factorial combinations of parameter inputs. 
### Read only
Checking "Read only" locks the *Scenario's* properties from being edited. 
### Merge Dependencies
Checking "Merge dependencies" allows the dependencies for the source *Scenarios* to be merged, rather than prioritized. For more information about dependencies, see [SyncroSim - Sharing Data](https://docs.syncrosim.com/how_to_guides/properties_dependencies.html)
### Ignore Dependencies...
Selecting "Ignore Dependencies" will lead to a pane showing Datafeeds and Packages. Checking the box "Ignore" allows the *Scenario* to ignore the selected package when running. 
<br>

<p id="heading02"> <h2>Pipeline</h2> </p>

The **Pipeline** datasheet controls the run order of the model transformers. 

### Stage
*Stages* represent the stages that will be run in the model. The general stages in **WISDM** include:
    1. Prepare Multiprocessing
    2. Spatial Data Preparation
    3. Data Preparation (Non-Spatial)
    4. Variable Reduction
    5a. Generalized Linear Model
    5b. Maxent
    5c. Random Forest
    6. Apply Model
Within this pane, the stages needed in the *Scenario* can be selected. 
### Run Order
*Run Order* represents the order in which the stages will be run. 
### Jobs
The *Jobs* argument is not currently used in WISDM, but within **SyncroSim**, it specifies the maximum number of jobs per stage. In **WISDM**, the number of jobs is defined by the number of multiprocessing tiles. 
<br>

<p id="heading03"> <h2>Datafeeds</h2> </p>

The **Datafeeds** datasheet controls which data are going into the *Scenario* and the properties of these data. 

### Data
The *Data* tab in this datasheet shows whethere the datafeed arguments have been filled in. When a new *Scenario* is opened, this tab will have nothing underneath it if no data or information in other *Datafeeds* have been input yet. If information <u>has</u> been input in these *Datafeeds*, a green checkmark will appear under this tab.
### Package Name
The *Package Name* tab specifies to which package each *Datafeed* belongs. For example, package names include *wisdm*, *corestime*, *core*, etc.
### View
The *View* tab will lead to the input datasheets for each *Datafeed*. Clicking the links under this tab will lead to the sections for the respective *Datafeeds*. 
### Datafeed
This tab shows the names of each *Datafeed* in the *Scenario*. The *Datafeeds* in WISDM include:
    Template Raster
    Covariate Data
    Field Data
    Field Data - Options
    Validation Options
    Site Data
    Covariate Selection Options
    Reduced Covariate List
    GLM
    Random Forest
    Maxent
    Model Outputs
    Output Options
### Source Scenario
The *Source Scenario* tab shows to which the *Scenario* the respecive *Datafeeds* and information belong. 
<br>
