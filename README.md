# Pseudogene

//  Created by Walfred MA in 2018, wangfei.ma@ucsf.edu.
//  Copyright © 2018 UCSF-Kwoklab. All rights reserved.

-----------description-----------

This is a pipeline to discovery and analyze all de novo retro-transcription events in de novo assemblies.

It has four pieces:

  1.find_retro.py	
  
  find all psuedo-like genes in de novo assemblies
  
  2.find_denovo.py	
  
  determine if psuedo-like genes are de novo or not by comparinfg with existing pseudogenes with anchors.
  
  3.locate_on_ref.py
  
  locate all psuedo-like genes on the reference genome
  
  4.analysis.py	
  
  analyze and determine if psuedo-like genes are retro-transcription events or not.
