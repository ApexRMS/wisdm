---
layout: default
section: 4
title: Variable Reduction
permalink: reference/variable-reduction
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
            <a href="#heading01"> &emsp;&emsp;&emsp;Site Data</a>
            <a href="#heading02"> &emsp;&emsp;&emsp;Covariate Selection Options</a>
            <a href="#heading03"> &emsp;&emsp;&emsp;Reduced Covariate List</a>
        {% endif %}
    {% endfor %}
</div>

# **Variable Reduction**

The **Variable Reduction** tab groups the following *Scenario Datasheets*::
* Site Data
* Covariate Selection Options
* Reduced Covariate List


In the SyncroSim UI, the **Variable Reduction** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

<br>

<p id="heading01"> <h2><b>Site Data</b></h2> </p>

The **Site Data** *Datasheet* contains information about covariate values at each *Field Data* site (*i.e.*, presence location or absence location provided in *Field Data*). 
> The *Site Data* *Datasheet* will be automatically populated after a *Scenario* has finished running. 

### **Site**
Sets the site IDs for each site containing data for a given covariate. The ID is based on the *Field Data* *Datasheet*.

### **Covariate**
Defines the covariate names for each site.

### **Value**
Sets the covariates values at each site.

<br>

<p id="heading02"> <h2><b>Covariate Selection Options</b></h2> </p>

The **Covariate Selection Options** *Datasheet* contains options for how correlations among covariates should be treated during the modeling process for *Stage 4 - Variable Reduction*.

### **Covariate Selection Method**
The *Covariate Selection Method* defines which method to use for variable reduction, and can either be Interactive (Correlation Viewer) or Automatic (Variance Inflation Factor). The Interactive method allows users to interact with the pop-up correlation viewer and select which variables to include in model development. The Automatic method automatically drops correlated variables using the "vifstep" function in the [usdm R package](https://cran.r-project.org/web/packages/usdm/usdm.pdf). 

<div class=indentation> 
    <i>Default:</i> Interactive (Correlation viewer).
</div>

### **Display Highest Correlations**
Determines whether only variables that are significantly correlated with the variable that has the <u>most total correlations</u> will be displayed in a correlation matrix. Other correlations between covariates, excluding the variable with the most total correlations, will not be displayed in the covariate correlation matrix. If "No" is selected, all correlations will be displayed in the matrix, including those not correlated with the variable that has the most correlations. 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Correlation Threshold**
Sets the threshold value for interpreting significant correlations between variables in the covariate correlation matrix interface. This argument should be a decimal between 0 and 1. For example, if the *Correlation Threshold* is set to 0.5, then any correlation cells with value 0.5 or higher will be colored. The gradient of correlation values will be colored from red to yellow, where red represents the highest correlation values and yellow represented the lower correlation values above the set threshold. 

<div class=indentation> 
    <i>Default:</i> 0.7.
</div>

### **Number of Plots**
The *Number of Plots* argument represents the number of covariates (an integer) that will be displayed in the covariate correlation matrix window. If "Yes" is selected for *Display Highest Correlations*, this argument should be left blank. This number can also be changed within the covariate correlation selection window while the *Scenario* is running. 

<div class=indentation> 
    <i>Default:</i> 5.
</div>

### **Variance Inflation Factor Threshold**
The *Variance Inflation Factor Threshold* is an optional argument that defines the variance inflation factor threshold for dropping correlated covariates, which is only used for the Automatic *Covariate Selection Method*.

<div class=indentation> 
    <i>Default:</i> 10.
</div>

<br>

<p id="heading03"> <h2><b>Retained Covariate List</b></h2> </p>

The **Retained Covariate List** *Datasheet* will be auto-populated once a *Scenario* has been run and covariates have been selected for modeling, and will show the covariates retained after manual or automatic covariate reduction. Users can optionally define variables here prior to variable reduction, and those variables will be prioritized and retained during automatic variable reduction.

### **Covariate**
The *Covariate* field shows which covariates have been included in the *Scenario*, and excludes covariates that were offered but removed from the *Scenario* during covariate selection.

<br>