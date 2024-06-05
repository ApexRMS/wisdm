---
layout: default
section: 5
title: Models
permalink: reference/models
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
            <a href="#heading01"> &emsp;&emsp;&emsp;GLM</a>
            <a href="#heading02"> &emsp;&emsp;&emsp;Random Forest</a>
            <a href="#heading03"> &emsp;&emsp;&emsp;Maxent</a>
            <a href="#heading04"> &emsp;&emsp;&emsp;Model Outputs</a>
        {% endif %}
    {% endfor %}
</div>

# **Models**

The **Models** tab groups the following *Scenario Datasheets*:
* Boosted Regression Tree
* GLM (Generalized Linear Model)
* Maxent
* Random Forest
* Model Outputs

In the SyncroSim UI, the **Models** tab can be accessed by right-clicking on a **WISDM** *Scenario* and selecting *Properties* from the context menu.

<br>

<p id="heading04"> <h2><b>Boosted Regression Trees</b></h2> </p>

The **Boosted Regression Trees** *Datasheet* contains information about the Boosted Regression Tree (BRT) algorithm, which is a machine-learning model and is described in and uses code from [Elith et al., 2008](https://doi.org/10.1111/j.1365-2656.2008.01390.x) and [Valavi et al., 2021](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1486). Arguments for this model correspond with model fitting inputs for the gbm.step function in the [dismo R package](https://cran.r-project.org/web/packages/dismo/dismo.pdf), and defaults were set based on Valavi et al., 2021.

### **Fitting Method**
The BRT fitting process can either be more user-specified or default-based. This argument allows the user to specify whether default values or user-specified values should be used in the fitting process for BRT. If this argument is set to "Use defaults and tuning", the code will tune the learning rate and number of trees to find values that fit a model to the full dataset and all CV splits.  If the argument is set to "Use values provided below", tuning does not occur and only the values provided in the below arguments will be used for model fitting.

<div class=indentation> 
    <i>Default:</i> Use defaults and tuning
</div>

### **Learning Rate**
A number between 0 and 1 that sets the weight applied to individual trees. A small learning rate means that each tree will have a smaller contribution to the overall model.

<div class=indentation> 
    <i>Default:</i> If not specified, learning rate will be determined based on the number of trees and tree complexity.
</div>

### **Number of trees added per stage**
Sets the number of initial trees, or trees added to the previously grown trees per round of model tuning.

<div class=indentation> 
    <i>Default:</i> 50.
</div>

### **Bag Fraction**
Sets the proportion of observations used in selecting variables.

<div class=indentation> 
    <i>Default:</i> 0.75.
</div>

### **Maximum number of trees**
Sets the maximum number of trees to fit before stopping.

<div class=indentation> 
    <i>Default:</i> 10,000.
</div>
<br>

<p id="heading01"> <h2><b>GLM</b></h2> </p>

The **GLM** *Datasheet* contains information about the Generalized Linear Model (GLM) algorithm options.

### **Select Best Predictors**
Determines whether to the GLM should use internal statistical methods to determine which predictors are significant to the model and  remove predictors that are unimportant for the *Scenario*. If *Select Best Predictors* is set to "No", the GLM will use all predictors offered to it. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Simplification Method**
Defines how the predictors are chosen for the GLM. Two options are available: "AIC" (Akaike Information Criterion) or "BIC" (Bayesian Information Criterion). AIC and BIC measure how well the GLM fits the data based on the covariates selected during the covariate selection step, and evaluates changes to the criteria when adding or dropping covariates. For more information about these criteria, see [Aho, Derryberry, & Peterson](https://doi.org/10.1890/13-1452.1).

<div class=indentation> 
    <i>Default:</i> AIC.
</div>

### **Consider Squared Terms**
Determines whether the GLM should consider relationships between the predictors and the response in a non-linear relationship by squaring the predictor values. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Consider Interactions**
Determines whether the GLM should take into account interactions among covariates during the modeling process. 

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
    <i>Default:</i> 10,000, or maximum number of raster cells.
</div>


### **Number of processing threads**
Sets the number of processing threads used for the Maxent algorithm. If this argument is not specified, the number of processing threads for Maxent will be set relative to the number of cores on the machine on which the **Scenario** is being run.

<div class=indentation> 
    <i>Default:</i> SyncroSim multiprocessing settings.
</div>

### **View maxent user interface during run**
Determine whether to show the Maxent user interface. 

<div class=indentation> 
    <i>Default:</i> No.
</div>


### **Save maxent formatted input/output files**
Determines whether to save Maxent's input and output files to the library folders.

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>

<p id="heading02"> <h2><b>Random Forest</b></h2> </p>

The **Random Forest** *Datasheet* contains information about the Random Forest (RF) algorithm - a machine-learning ensemble classifier. More information about Random Forest can be found in the [randomForest R package documentation](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf).

### **Evaluate covariate importance**
Determines whether the RF algorithm should evaluate covariate importance during the modeling process by calculating the change in fit statistics for trees that include each covariate. 

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Calculate casewise importance**
Determines whether the RF algorithm should calculate the importance of each covariate during the classification process and measure how significantly covariates influence the output. For more information about covariate and casewise importance in the RF algorithm, see [Classification and Regression with Random Forest](https://haoen-cui.github.io/SOA-Exam-PA-R-Package-Documentation/randomForest/reference/randomForest.html). 

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
Determines whether each decision tree's vote should be divided by the total number of votes to ensure that the sum of the votes is 1. Votes will be expressed as fractions.

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

### **Calculate proximity**
Determines whether proximity measures should be calculated. Proximity measures are based on the frequency that pairs of data points are in the same terminal nodes.

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Sample with replacement**
Determines whether each decision tree should be trained on a bootstrap sample of the data, so some observations can be repeated or left out.

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br>


<!-- <p id="heading05"> <h2><b>Generalized Additive Model</b></h2> </p>

The **Generalized Additive Model** is an algorithm similar to the GLM model but allows for nonlinearity in the fitted functions and is described in [Valavi et al., 2021](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1486). 

### **Allow smooth terms to shrink to zero**
If set to "Yes", sets smooth terms to shrinkage smoothers in which a small multiple of the identity matrix is added to the smoothing penalty, so that strong enough penalization will shrink all coefficients of the smooth to zero. This approach allows smoother terms to effectively be penalized out of the model during smoothing parameter estimation. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Consider linear terms**
If set to "Yes", linear terms for each variable are added to the model in addition to the smooth terms. 

<div class=indentation> 
    <i>Default:</i> No.
</div>

<br> -->

<p id="heading06"> <h2><b>Model Outputs</b></h2> </p>

The **Model Outputs** *Datasheet* contains information about the outputs of the models. These outputs will be visible in the **Maps** tab in the lower left panel when the **Scenario's** results have been added and can be exported from the **Export** tab. The outputs may include: 
* Response Curves
* Standard Residuals Plots
* Residual Smooth Plot
* Residual Smooth RDS
* Text Output
* Calibration Plot
* ROC/AUC Plot
* AUCPR Plot
* Confusion Matrix
* Variable Importance Plot
* Variable Importance Data
* Maxent Files

Note that not all outputs are generated per model run. The outputs that are generated depend on the type of model selected and the type of field data being used (e.g., presence/absence, count).

### **Model RDS**
Defines the names of the model .rds files (which hold R objects) for algorithms that have completed execution.

### **Response Curves**
Model **Response Curves** show the relationship between each predictor included in the model and the fitted values of the model. In general, these response curves should be smooth positive and negative relationships in agreement with the biological relationships that occur in the environment. While some exceptions occur, bumpy response curves can indicate the need for an algorithm to be optimized.

### **Standard Residuals Plots**
Model **Residuals plots** show the plotted model deviance residuals for the algorithm.

### **Residual Smooth Plot**
Model **Residual Smooth Plots** show the spatial relationship between the model deviance residuals. With the assumption that the residuals will be independent from each other, a spatial pattern in the deviance residuals could indicate an issue with the model fit, and these patterns can be seen in clusters of higher or lower residuals [(Dormann et al., 2007)](https://doi.org/10.1111/j.2007.0906-7590.05171.x).

### **Residual Smooth RDS**
Contains a file path to the .rds file holding the information pertaining to the **Residual Smooth Plots**.

### **Text Output**
The **Text Outputs** contain information about the algorithm's settings and results. For example, the text output for **GLM** contains information about the model family and simplification method used as well as the number of covariates in the final model and evaluation statistics, along with other information.

### **Calibration Plot**
The **Calibration Plot** shows the relationship between the predicted values from the model and the actual observations from the test proportion of the data. Calibration plots showing a higher agreement between the predicted probability of presence and the probability of presence indicate good calibration of the model, but the AUC of the model should be considered along with the calibration plot [(Pearce & Ferrier, 2000;](https://doi.org/10.1016/S0304-3800(00)00322-7)[ Elith et al., 2010)](https://doi.org/10.1111/j.2041-210X.2010.00036.x).

### **ROC/AUC Plot**
The **ROC/AUC Plots** show the relationship between sensitivity (true positives) and specificity (False positives) in the model algorithm. The output of these curves depend on whether the *Split Data For Model Training and Testing* and *Use Cross Validation for Model Selection* arguments were set to "Yes" in the **Validation Options** datasheet in the **Data Preparation** tab. It shows the sensitivity versus specificity of the training split and testing split or the training split and cross validation mean, along with AUC values. Better curves are generally those that arch far above the diagonal of the plot and those that have smaller discrepancies between the training and testing/validation data. 

### **AUCPR Plot**
The **AUCPR** plots show the relationship between recall (True positives / (True positives + False negatives)) and precision (True positives / (True positives + False Positives)) of the algorithm along with the Training split AUC and cross validation mean AUC. For more information about ROC/AUC and AUCPR plots, see [(Davis & Goadrich, 2006)](https://www.researchgate.net/publication/215721831_The_Relationship_Between_Precision-Recall_and_ROC_Curves#full-text). 

### **Confusion Matrix**
The **Confusion Matrix** Shows the number of values classified by the algorithm as a presence or an absence compared to the observed number of presences and absences in the data. The algorithm will output a confusion matrix for the train data, and a confusion matrix for the cross validation or test data. Each matrix consists of 4 cells, which predicted values on the left and observed values on the bottom. The top left cell shows the number of values predicted to be a presence that were observed to be a presence. The top right cell shows the number of values predicted to be a presence that were actually observed as absences. The bottom left cell shows the number of values predicted to be an absence that were actually observed as presences. The bottom right cell shows the number of values predicted to be an absence that were also observed to be absences. These values contribute to the calculation of the statistics at the bottom of the matrix: percent correctly classified, sensitivity, specificity, true skill statistic, and Cohen's Kappa. The values in the **Confusion Matrix** can help identify how well or poorly the algorithm has made its predictions.

### **Variable Importance Plot**
The **Variable Importance Plot** shows the relative influence of each predictor in the model. Variable importance is calculated by permuting the values of the predictors in the dataset and predicting the model to the new dataset with the permuted predictor values and measuring the mean drop in the AUC value using 5 different permutations of the predictor. The importance is measured by the change in AUC when each predictor is permuted and appears on the x-axis with individual variables on the y-axis. The importance is shown for the cross-validation, train, and test data.

### **Variable Importance Data**
Exporting the **Variable Importance Data** outputs a .csv file containing variable importance values for the train data, test data, and each cross validation split.

### **Maxent Files** 
The **Maxent Files** are files that are output while the Maxent algorithm is running and are output in a .zip folder. 

<br>