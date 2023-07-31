---
layout: default
section: 04
title: Models
permalink: section04
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

This is the page for the **Models** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

<br>
<br>

<p id="heading01"> <h2>GLM</h2> </p>

The **GLM** datasheet contains information about the Generalized Linear Model algorithm options.

### Select Best Predictors
Selecting "Yes" in the *Select Best Predictors* argument means that the GLM will use internal statistical methods to determine which predictors are significant to the model and will remove predictors that are unimportant for the *Scenario*. Selecting "No" requires the GLM to use all predictors offered to it.
### Simplification Method
The *Simplification Method* will either be AIC (Akaike Information Criterion) or BIC (Bayesian Information Criterion), and these criteria relate to how the predictors are chosen for the GLM. For more information about these criteria, see [Aho, Derryberry, & Peterson](https://doi.org/10.1890/13-1452.1).
### Consider Squared Terms
Selecting "Yes" for the *Consider Squared Terms* argument allows the GLM to consider relationships between the predictors and the response in a non-linear relationship by squaring the predictor values. 
### Consider Interactions
Selecting "Yes" for the *Consider Interactions* argument will allow the GLM to take interactions among covariates into account during the modeling process.
<br>

<p id="heading02"> <h2>Random Forest</h2> </p>

The **Random Forest** datasheet contains information about the Random Forest (RF) algorithm options. More information about random forest can be found in the [CRAN - Random Forest Reference Manual](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf).

### Evaluate covariate importance
Selecting "Yes" for the *Evaluate covariate importance* argument will allow the RF algorithm to evaluate covariate importance during the modeling process. 
### Calculate casewise importance
Selecting "Yes" for the *Calculate casewise importance* argument will allow RF to calculate the importance of each covariate during the classification process, and measures how significantly covariates influence the output. For more information about covariate and casewise importance in the RF algorithm, see [Classification and Regression with Random Forest](https://haoen-cui.github.io/SOA-Exam-PA-R-Package-Documentation/randomForest/reference/randomForest.html)
### Number of variables sampled at split
The *Number of variables sampled at split* argument will auto-populate after the model runs, and shows the number of variables actually used in the random forest, if not specified beforehand.
### Maximum number of nodes
The *Maximum number of nodes* argument specifies the maximum number of nodes the RF algorithm can have.
### Number of trees
If not specified, the *Number of trees* argument will auto-populate after the RF algorithm has finished running. This argument controls the number of trees in the random forest, and the default is 1000.
### Node size
The *Node size* argument specifies the size of the terminal nodes. A larger node size causes smaller trees to be grown and will take less time. This argument can also auto-populate after the RF algorithm has finished running, with a default node size of 1.
### Normalize votes
Selecting "Yes" on the *Normalize votes* argument means that each decision tree's vote is divided by the total number of votes to ensure that the sum of the votes is 1 (votes will be expressed as fractions). The default for this argument is "Yes".
### Calculate proximity
Selecting "Yes" on the *Calculate proximity* argument means that proximity measures based on the frequency that pairs of data points are in the same terminal nodes will be calculated. The default for this argument is "No".
### Sample with replacement
Selecting "Yes" on *Sample with replacement* will allow each decision tree to be trained on a bootstrap sample of the data, so some observations can be repeated or left out. The default for this argument is "No".

<br>
<p id="heading03"> <h2>Maxent</h2> </p>

The **Maxent** datasheet contains information about the Maxent algorithm options.

### Memory allocation (GB)
The *Memory allocation (GB)* argument represents the memory that should be allocated to the Maxent algorithm, in gigabytes. If left blank, this field will auto-populate after the algorithm is finished running.
### Maximum number of background points
The *Maximum number of background points* argument specifies the maximum number of background points that will be generated during the modeling process for the Maxent algorithm. If left blank, this field will auto-populate at the end of the modeling process. More background points will lead to a higher entropy of the model, but only one background point can be generated per raster cell. The default for this argument is 10,000. 
### Number of processing threads
The *Number of processing threads* argument specifies the number of processing threads used for the Maxent algorithm. If this argument is not specified, the number of processing threads for Maxent will rely on the multiprocessing environment settings in **SyncroSim**.
### View maxent user interface during run
Selecting "Yes" on the *View maxent user interface during run* argument will show the Maxent user interface. The default for this argument is "No", as showing the user interface will often take more processing time and power. 
### Save maxent formated input/output files
Selecting "Yes" on the *View maxent formatted input/output files* argument will allow Maxent input and output files to be saved to the library folders. 
<br>

<p id="heading04"> <h2>Model Outputs</h2> </p>

The **Model Outputs** datasheet contains information about the outputs of the models.

### Model RDS
This table will populate after models have been run, and show the names of the model .rds files for algorithms that have completed execution.
<br>

