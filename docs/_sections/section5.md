---
layout: default
section: 5
title: Models
permalink: reference/models
---

<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;GLM</a>
            <a href="#heading02"> &emsp;Random Forest</a>
            <a href="#heading03"> &emsp;Maxent</a>
            <a href="#heading04"> &emsp;Model Outputs</a>
        {% endif %}
    {% endfor %}
</div>

# **Models**

The **Models** tab groups the following *Scenario Datasheets*:
* GLM
* Random Forest
* Maxent
* Model Outputs

In the SyncroSim UI, the **Models** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

<br>

<p id="heading01"> <h2><b>GLM</b></h2> </p>

The **GLM** *Datasheet* contains information about the Generalized Linear Model (GLM) algorithm options.

### **Select Best Predictors**
Determines whether to the GLM should use internal statistical methods to determine which predictors are significant to the model and  remove predictors that are unimportant for the *Scenario* ("Yes"). If *Select Best Predictors* is set to "No", the GLM will use all predictors offered to it. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Simplification Method**
Defines how the predictors are chosen for the GLM. Two options are available: "AIC" (Akaike Information Criterion) or "BIC" (Bayesian Information Criterion). AIC and BIC measure how well the GLM fits the data based on the covariates selected during the covariate selection step, and evaluates changes to the criteria when adding or dropping covariates. For more information about these criteria, see [Aho, Derryberry, & Peterson](https://doi.org/10.1890/13-1452.1).

<div class=indentation> 
    <i>Default:</i> AIC.
</div>

### **Consider Squared Terms**
Determines whether the GLM should consider relationships between the predictors and the response in a non-linear relationship by squaring the predictor values ("Yes"). 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Consider Interactions**
Determines whether the GLM should take into account interactions among covariates during the modeling process ("Yes"). 

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>

<p id="heading02"> <h2><b>Random Forest</b></h2> </p>

The **Random Forest** *Datasheet* contains information about the Random Forest (RF) algorithm - a machine-learning ensemble classifier. More information about Random Forest can be found in the [randomForest R package documentation](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf).

### **Evaluate covariate importance**
Determines whether the RF algorithm should evaluate covariate importance during the modeling process by calculating the change in fit statistics for trees that include each covariate ("Yes"). 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Calculate casewise importance**
Determines whether the RF algorithm should calculate the importance of each covariate during the classification process and measure how significantly covariates influence the output ("Yes"). For more information about covariate and casewise importance in the RF algorithm, see [Classification and Regression with Random Forest](https://haoen-cui.github.io/SOA-Exam-PA-R-Package-Documentation/randomForest/reference/randomForest.html). 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Number of variables sampled at split**
*Optional*. Sets the number of variables randomly sampled as candidates in each split of the random forest.

### **Maximum number of nodes**
*Optional*. Sets the maximum number of nodes the RF algorithm can have.

### **Number of trees**
Sets the number of trees for the RF algorithm. 

<div class=indentation> 
    <i>Default:</i> 1000.
</div>

### **Node size**
*Optional*. Sets the size of the terminal nodes. A larger node size causes smaller trees to be grown and will take less time.

### **Normalize votes**
Determines whether each decision tree's vote should be divided by the total number of votes to ensure that the sum of the votes is 1 ("Yes"). Votes will be expressed as fractions.

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Calculate proximity**
Determines whether proximity measures should be calculated ("Yes"). Proximity measures are based on the frequency that pairs of data points are in the same terminal nodes.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Sample with replacement**
Determines whether each decision tree should be trained on a bootstrap sample of the data, so some observations can be repeated or left out ("Yes").

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>

<p id="heading03"> <h2><b>Maxent</b></h2> </p>

The **Maxent** *Datasheet* contains information about the Maxent algorithm, a machine learning technique that minimizes relative entropy between the probability density of the species and of the environment.

### **Memory allocation (GB)**
Sets the memory that should be allocated to the Maxent algorithm, in gigabytes.

<div class=indentation> 
    <i>Default:</i> 512 GB.
</div>

### **Maximum number of background points**
Sets the maximum number of background points that will be generated during the modeling process for the Maxent algorithm. More background points will lead to a higher entropy of the model, and only one background point can be generated per raster cell. 

<div class=indentation> 
    <i>Default:</i> 10,000.
</div>


### **Number of processing threads**
Sets the number of processing threads used for the Maxent algorithm. If this argument is not specified, the number of processing threads for Maxent will rely on the .

<div class=indentation> 
    <i>Default:</i> SyncroSim multiprocessing settings.
</div>

### **View maxent user interface during run**
Determine whether to show the Maxent user interface ("Yes"). 

<div class=indentation> 
    <i>Default:</i> No.
</div>


### **Save maxent formated input/output files**
Determines whether to save Maxent's input and output files to the library folders ("Yes").

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>

<p id="heading04"> <h2><b>Model Outputs</b></h2> </p>

The **Model Outputs** *Datasheet* contains information about the outputs of the models.

### **Model RDS**
Defines the names of the model .rds files (which hold R objects) for algorithms that have completed execution.

<br>