---
layout: default
title: Reference
description: "Reference guide for the WISDM package"
permalink: /reference
---

# Reference guide for the **WISDM** SyncroSim *Package*

In SyncroSim, all of the of inputs and outputs associated with a model are stored in a single file (with the extension .ssim), referred to as a SyncroSim *Library*. The inputs and outputs contained within a SyncroSim *Library* are organized into *Datafeeds* (or tabs in the SyncroSim UI). Each *Datafeed* can be made up of one or more tables of data, called *Datasheets*. Each *Datasheet* is associated with one of three scopes: *Library*, *Project* or *Scenario*. See the [SyncroSim documentation](https://docs.syncrosim.com/how_to_guides/library_overview.html) for more details.

**WISDM** has a wide variety of *Project* and *Scenario* scoped *Datasheets*. Within this Reference guide, you will find details on the function of each field within **WISDM** *Datasheets*.

**WISDM** *Datasheets* are organized under the following tabs:

{% for section in site.sections %}
  <li> <a href="{{site.baseurl}}{{ section.url }}"> {{ section.title }}</a> </li>
{% endfor %}

<br>