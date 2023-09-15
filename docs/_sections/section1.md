---
layout: default
section: 1
title: Covariates
permalink: reference/covariates
---

<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
        {% endif %}
    {% endfor %}
</div>

# **Covariates**

The **Covariates** tab is a *Project Datasheet*.

In the SyncroSim UI, it can be accessed by right-clicking on a **WISDM** *Project*, selecting *Properties* from the context menu, and navigating to the *Covariates* tab.

<br>

### **Covariate Name**
Provides a name for the covariate. The *Covariate Name* is used by the **Site Data** *Datasheet* under the [**Variable Reduction**](variable-reduction#heading01) tab.  

### **Is Categorical**
*Optional*. Determines whether to the covariate is factor-based ("Yes"). 

<div class=indentation> 
    <i>Default:</i> No.
</div>

### **Description**
*Optional*. Provides a description of the covariate.

### **ID**
*Optional*. Sets the site IDs for each covariate.

### **Color**
*Optional*. Defines a color for each covariate.

<br>