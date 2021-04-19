
.. _intro:

Introduction to roboscheduler
===============================

``roboscheduler`` is the scheduling utility for SDSS V FPS program. The goal of ``roboscheduler`` is to schedule the most "optimal" or "important" field at any given time, based on how difficult it is to observe that field (e.g. if it needs 8 exposures vs 1) and cadence requirements.

Roboscheduler context
---------------------

During the FPS program, the ``roboscheduler`` product will be used almost exclusively with the web app scheduling tool ``Kronos`` that is currently in development. Normal users will rarely need to interact directly with the scheduler. 

In order to optimize the scheduling strategy, simulations using the scheduler are regularly performed using the ``observesim`` package. This package mimics stepping through time and interacting with the scheduler much the same way ``Kronos`` will.
