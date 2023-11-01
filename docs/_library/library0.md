---
layout: default
section: 0
title: Network
permalink: reference/network
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
        {% endif %}
    {% endfor %}
</div>

# **Network**

The **Network** node is a *Library Datasheet*.

In the SyncroSim UI, it can be accessed by right-clicking on a **WISDM** *Library*, selecting *Properties* from the context menu, and navigating to the *Options* tab.

<br>

### **Access online transformation library**
Determines whether to use [PROJ](https://proj.org/en/9.3/) to query its online repository for the best transformation available for projecting spatial data between two different CRS ("Yes"). If "No" is selected, then PROJ is only able to use transformations stored in the **WISDM** *Library*, which get downloaded when PROJ is installed. Enabling network access is recommended, as depending on the desired transformation, a more accurate option maybe available online. Nevertheless, disabling network access is an option for handling restricted network issues or for those who wish to run **WISDM** offline.

> [PROJ](https://proj.org/en/9.3/) is used to transform spatial data in the "Spatial Data Preparation" [*Stage*](general#heading02).  

<div class=indentation> 
    <i>Default:</i> Yes.
</div>

<br>