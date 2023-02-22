---
layout: default
title: Reference
description: "Reference guide for the WISDM package"
permalink: /reference
---

# Reference for the **WISDM** SyncroSim Package

A SyncroSim model can have any number of Properties, each of which can have a *scope* of **Library**, **Project**, or **Scenario**. **WISDM** has a wide variety of Project and Scenario scoped Properties organized into *datasheets*. *Datasheet*, in turn, are organised under various tabs in the **SyncroSim** user interface. 

Within this Reference guide, you will find details on the function of each column within each datasheet. *Datasheets* are organized under the following tabs:

{% for section in site.sections %}
  <li> <a href="{{site.baseurl}}{{ section.url }}"> {{ section.title }}</a> </li>
{% endfor %}