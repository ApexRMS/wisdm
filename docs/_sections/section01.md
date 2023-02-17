---
layout: default
section: 01
title: General
permalink: section01
---


<!--- Sidebar Navigation Menu --->
<div class="sidenav">
    {% for section in site.sections %}
        <a href="{{ section.url }}"> <b>{{ section.title }}</b> </a>
        {% if section.section == page.section %}
            <a href="#heading01"> &emsp;Summary</a>
            <a href="#heading02"> &emsp;Pipeline</a>
            <a href="#heading03"> &emsp;Datafeeds</a>
        {% endif %}
    {% endfor %}
</div>
<br>

# **General**

This is the page for the **General** tab with content organized into headings and, optionally, subsections.

The **Sidebar Navigation Menu** lists all the headings within the page, and all the Sections within the Reference guide. 
<br>
<br>

<p id="heading01"> <h2>Summary</h2> </p>

The **Summary** datasheet contains general information about the *Scenario*. 

### Name
This column details the name of the *Scenario*.
### Owner

### Description

### Project

### Library

### Last modified

### Auto generation tags

### Read only

### Merge Dependencies

### Ignore Dependencies...

<br>

<p id="heading02"> <h2>Pipeline</h2> </p>

The **Pipeline** datasheet controls the run order of the model transformers. 

### Stage

### Run Order

### Jobs

<br>

<p id="heading03"> <h2>Datafeeds</h2> </p>
Some random text here just to show the structure of the page and how internal links work.

### Data

### Package Name

### View

### Datafeed

### Source Scenario

<br>
