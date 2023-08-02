---
layout: default
section: 03
title: Variable Reduction
permalink: section03
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
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

The **Site Data** datasheet contains information about covariate values at each *Field Data* site (i.e. presence location or absence location provided in *Field Data*). This datasheet will auto-populate after a *Scenario* has finished running. 

### Site
The *Site* column contains site IDs from the *Field Data* datasheet within the **Data Preparation** tab. It is repeated within this column for each covariate containing data for this site.
### Covariate
The *Covariate* column contains the covariate names for each site.
### Value
The *Value* tab contains the values of covariates at each site. 
<br>

<p id="heading02"> <h2>Covariate Selection Options</h2> </p>

The **Covariate Selection Options** datasheet contains options for how correlations among covariates should be treated during the modeling process for *Stage 4 - Variable Reduction*.

### Display Highest Correlations
Selecting "Yes" for *Display Highest Correlations* means that only the variables that are significantly correlated with the variable that has the <u>most total correlations</u> will be displayed in a correlation matrix. Other correlations between covariates, excluding the variable with the most total correlations, will not be shown in the covariate correlation matrix. If "No" is selected, this means that all correlations will be shown in the matrix, including those not correlated with the variable that has the most correlations.
### Correlation Threshold
The *Correlation Threshold* argument relates to the covariate correlation matrix interface. This argument should be a decimal between 0 and 1. For example, if the *Correlation Threshold* is 0.5, then any correlation cell value at 0.5 or higher will be colored to draw attention to the correlation. Higher correlations will be represented with redder shades, and lower correlations that are still above the threshold will be represented with lighter yellow shades. 
### Number of Plots
The *Number of Plots* argument represents the number of covariates (an integer) that will be shown in the covariate correlation matrix window. If "Yes" is selected for *Display Highest Correlations*, this argument should be left blank. This number can also be changed within the covariate correlation selection window while the *Scenario* is running. 

<br>

<p id="heading03"> <h2>Reduced Covariate List</h2> </p>

The **Reduced Covariate List** will be auto-populated once a *Scenario* has been run and covariates have been selected for modeling. 

### Covariate
The *Covariate* field shows which covariates have been included in the *Scenario*, and excludes covariates that were offered but removed from the *Scenario* during covariate selection.
<br>

