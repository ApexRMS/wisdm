---
layout: default
section: 6
title: Output Options
permalink: reference/output-options
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
            <a class="indent1" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% else %}
            <a class="selected" href="{{site.baseurl}}{{ section.url }}"> &emsp;{{ section.title }} </a>
        {% endif %}
    {% endfor %}
</div>

# **Output Options**

The **Output Options** tab is a *Scenario Datasheet*.

In the SyncroSim UI, the **Output Options** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

For more information about viewing outputs, see the [SyncroSim documentation](https://docs.syncrosim.com/how_to_guides/results_overview.html). Information about specific **Map Outputs**, **Model Outputs** and **Ensemble Outputs** can also be found below.

<br>

## **Output Options**

### **Make Probability Map**
Determines whether to create a probability map with the extent of the template raster for the **Scenario**. 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Make Residuals Map**
Determines whether to create a model deviance residuals map for the **Scenario**.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make Multivariate Environmental Similarity Surface (MESS) Map**
Determines whether to create a MESS map for the **Scenario**.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make Most Dissimilar Variable (MoD) Map**
Determines whether to create a MoD map for the **Scenario**.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make Binary Map**
Determines whether to create a binary map, which is generated from the probability map using thresholds calculated from the model.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Binary Threshold Optimization Method**
The **Binary Threshold Optimization Method** argument indicates which method do use for setting the threshold value when generating the binary map. Options include:
* <u>Max kappa</u>: The threshold at which kappa is highest.
* <u>Max sensitivity and specificity</u>: The threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest.
* <u>No omission</u>: The highest threshold at which there is no omission.
* <u>Prevalence</u>: The threshold at which modeled prevalence is closest to observed prevalence.
* <u>Sensitivity equals specificity</u>: The threshold at which there is equal sensitivity and specificity.

<br>

## **Map Outputs**

The **Map Outputs** datasheet contains information about the spatial outputs of the **Scenario**. The spatial outputs can be found in the **Maps** tab in the lower left panel when the **Scenario's** results have been added and can be exported from the **Export** tab. The outputs may include:
* Probability Map
* MESS Map
* MoD Map
* Residuals Map

Note that not all map outputs are generated per model run. Whether they are generated is determined by the arguments in the **Output Options** datasheet. 

### **Probability Map**
The **Probability Map** is the main output of the fitted model and shows probability values based on input **Field Data** and **Covariate Data**. These probability values represent the percent likelihood of occurrence in different areas of the **Template Raster** and generally range from values of 0 - 100%.

### **MESS Map**
The **MESS Map** is the Multivariate Environmental Similarity Surface, which represents values as positive, negative, or zero. This map shows how well locations on the **Template Raster** fit into the range of covariate data to which the training data were fit. Positive areas on this map represent areas where the covariate ranges are more similar to those to which the training data of the model were fit. Negative areas on this map represent areas where the covariate ranges are not similar to those to which the training data of the model were fit. Values of zero on this map represent areas where ranges of covariate data at these locations and ranges of covariate data to which the training data of the model were fit are marginally similar [(Elith et al., 2010)](https://doi.org/10.1111/j.2041-210X.2010.00036.x).

### **MoD Map**
The **MoD Map** is a map of the most dissimilar variable. This map is similar to the **MESS Map** in that it shows regions where covariate ranges were most dissimilar from those used to fit the training data. However, this map shows which covariates used in the model was furthest from the range of the observations used for model training and where.

### **Residuals Map**
The **Residuals Map** is a similar output to the **Residuals Smooth Plot** in the **Model Outputs** datafeed, but is mapped to the extent of the **Template Raster** rather than on an x/y plot of coordinates. This map shows the spatial relationship between the model deviance residuals. A spatial pattern in the model deviance residuals could indicate an issue with the model fit, and can be identified through spatial clusters of high or low residuals.

<br>

## **Ensemble Options**
The **Ensemble Options** datasheet defines if and how the models should be ensembled to integrate outputs from multiple algorithms. 

### **Make probability ensemble map**
The **Make probability ensemble map** argument defines whether an ensemble map should be created for the scenario. This would combine provided probability maps into a single probability ensemble map. The ensemble creation process is impacted by the provided inputs for the **Probability ensemble method**, **Normalize probability before ensemble**, and **Ignore predicted NA values** arguments. 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Probability ensemble method**
The **Probability ensemble method** argument determines the method for ensembling and can either be a mean or a sum of probability values from different algorithms. 

<div class=indentation> 
    <i>Default:</i> Mean
</div>

### **Normalize probability before ensemble**
The **Normalize probability before ensemble** argument specifies whether each probability map should be normalized before ensembling. If set to "Yes", each probability map will be normalized so that values range between 0 and 100 before combining the map information for ensemble output. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Make binary ensemble map**
If the **Make binary ensemble map** argument is set to "Yes", binary maps will be combined to create a single binary ensemble map. Provided inputs for the **Binary ensemble method** and **Ignore predicted NA values** aruments will further impact the binary ensemble map created. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Binary ensemble method**
The **Binary ensemble method** argument defines the ensembling method as either "Mean" or "Sum", so the ensemble will be created using a mean or sum of the inputs. The mean or sum of multiple binary inputs will produce an ensemble that is not binary. 

 <div class=indentation> 
    <i>Default:</i> No.
</div>

### **Ignore predicted NA values**
If the **Ignore predicted NA values** argument is set to "Yes", NA values will be ignored during ensemble calculation (i.e., if one input map has an NA value for a given cell, this value will be ignored and the output value for this cell is calculated from the cell's value in other input maps). If the **Ignore predicted NA values** argument is set to "No", NA values will be included in ensemble calculation, meaning that if one input map has an NA value in a given cell, the output ensemble map will have an NA value for this cell. 

 <div class=indentation> 
    <i>Default:</i> Yes.
</div>

<br>

## **Ensemble Outputs**

The **Ensemble Outputs** datasheet contains information about the ensembled outputs from the **Scenario**. The ensemble outputs can be found in the **Maps** tab in the lower left panel when the **Scenario's** results have been added and can be exported from the **Export** tab. The outputs may include:
* Probability Ensemble (Mean)
* Probability Ensemble (Sum)
* Binary Ensemble (Mean)
* Binary Ensemble (Sum)

Note that only one probability and one binary ensemble can created for the result **Scenario**.  The ensemble output type will depend on the inputs (or defaults) in the **Ensemble Outputs** datasheet.

### **Probability Ensemble (Mean)**
The **Probability Ensemble (Mean)** is a map of probability created by taking the per-pixel average of values in the **Scenario**'s individual probability outputs. Values represent the likelihood of occurrence ranging from 0 - 100%.

### **Probability Ensemble (Sum)**
The **Probability Ensemble (Sum)** is a map created by taking the per-pixel sum of values in the **Scenario**'s individual probability outputs.

### **Binary Ensemble (Mean)**
The **Binary Ensemble (Mean)** is a map of probability created by taking the per-pixel average of values in the **Scenario**'s individual binary outputs.

### **Binary Ensemble (Sum)**
The **Binary Ensemble (Sum)** is a map created by taking the per-pixel sum of values in the **Scenario**'s individual binary outputs.  
