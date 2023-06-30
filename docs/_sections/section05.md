---
layout: default
section: 05
title: Output Options
permalink: section05
---

<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{site.baseurl}}{{ section.url }}"> <b>{{ section.title }}</b> </a>
    {% endfor %}
</div>

# **Output Options**

This is the page for the **Output Options** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 

### Make Probability Map
Selecting "Yes" on the *Make Probability Map* argument will create a probability map for the **Scenario**.
### Make Residuals Map
Selecting "Yes" on the *Make Residuals Map* argument will create a residuals map for the **Scenario**.
### Make Multivariate Environmental Similarity Surface (MESS) Map
Selecting "Yes" on the *Make Multivariate Environmental Similarity Surface (MESS) Map* argument will create a MESS map for the **Scenario**.
### Make Most Dissimilar Variable (MoD) Map
Selecting "Yes" on the *Make Most Dissimilar Variable (MoD) Map* argument will create a MoD map for the **Scenario**.