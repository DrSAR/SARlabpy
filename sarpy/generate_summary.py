# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:10:05 2013

@author: stefan
"""
import re
import sarpy
import os
import collections
from collections import defaultdict
import json
import numpy

def generate_summary(masterlist_name):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)
        
    # Generate the filename of the output file based on the masterlist
    filename = root + '_summary.pdf'

    ## Construct the data for the checks table
    
    setvals = set()
    pats = []
    for k,v in master_list.items():
        pats.append(k)
        for val in v:
            setvals.add(val)
    
    # The order of this setvals is a bit weird, so reverse the set, and turn it into a list:
    headings = natural_sort(list(setvals))
    headings.insert(0,'')
    
    checks =[]
    checks.append(headings)
    for p in pats:
        curr=[]
        curr.append(p)
    
        #Recall the first one is blank to match up with the patients column
        for k in headings[1:]:
            
            if master_list[p][k][0] == '':
                curr.append('No')
            else:
                curr.append('Yes')
        checks.append(curr)
    
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.pagesizes import letter, inch, landscape, cm 
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, PageBreak
    from reportlab.pdfgen import canvas


    doc = SimpleDocTemplate(filename, pagesize=landscape(letter))

    width, height = 8.5,14
    # container for the 'Flowable' objects
    elements = []
     
    checks_size = numpy.array(checks).shape

    #return checks
    colwidths = checks_size[1]*[0.7*inch]
    colwidths.insert(0,0.95*inch)
    
    rowHeights = checks_size[0]*[0.3*inch]
    
    t=Table(checks,colwidths, rowHeights)
    styled = []
    styled = [ ('TEXTCOLOR',(0,0),(0,-1),colors.blue),
               ('TEXTCOLOR',(0,0),(-1,0),colors.blue),
               ('ALIGN',(0,0),(-1,-1),'CENTER'),
               ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
               ('INNERGRID', (1,1), (-1,-1), 0.25, colors.black),
               ('BOX', (1,1), (-1,-1), 0.25, colors.black),
               ('BACKGROUND',(1,1),(-1,-1),colors.lightgreen)
            ]
    
    for cols in range(checks_size[0]):
    
        indexes = [i for i,x in enumerate(checks[cols]) if x == 'No']
        
        for idx in indexes:
            styled.append(('BACKGROUND',(idx,cols),(idx,cols),colors.orange))
            
    
    t.setStyle(TableStyle(styled))
    elements.append(t)

    elements.append(PageBreak())

    ## Construct the table for the scan information

    combinedParams = []

    styles = getSampleStyleSheet()

    for stype in natural_sort(list(setvals)):

        scanType = set()
        fov = set()
        px_size = set()
        scanParams = set()
        adata = set()

        for k,v in master_list.items():

            if master_list[k][stype][0] != "":
                scn = sarpy.Scan(master_list[k][stype][0])
                pack_extent = [scn.method.PVM_SPackArrSliceDistance[0]*scn.method.PVM_SPackArrNSlices[0] - scn.method.PVM_SPackArrSliceGap[0]]

                scanType.add(scn.method.Method)
                fov.add(str(scn.method.PVM_FovCm+pack_extent).strip('[,]'))
                px_size.add(str(scn.method.PVM_SpatResol+[scn.method.PVM_SliceThick]).strip('[,]'))
                scanParams.add(str([scn.method.PVM_RepetitionTime,scn.method.PVM_EchoTime1]).strip('[,]'))
                adata.add(str(list(scn.adata.keys())).strip('[,]'))

        if len(fov) != 1:
            fov = [str(len(fov)) + ' diff']

        if len(px_size) != 1:
            px_size = [str(len(px_size)) + ' diff']

        if len(scanParams) != 1:
            scanParams = [str(len(scanParams)) + ' diff']

        combinedParams.append([stype, list(fov)[0],list(px_size)[0],list(scanParams)[0],Paragraph(list(adata)[0],styles['BodyText'])])
    
    combinedParams.insert(0,['','FoV (x,y,z) [mm]','SpatRes (x,y,z) [mm]','TR, TE [ms]', 'Available adata'])  

    ##This transposes the list
    #combinedParams = map(None,*combinedParams) 

    #Ignore the list-->array-->list conversion for now... it's so that I can transpose things safely
    # Yes, there's probably a better way, go figure it out if you're so clever
         
    checks_size = numpy.array(combinedParams).shape

    #return combinedParams

    colwidths = (checks_size[1]-1)*[1.4*inch]
    colwidths.append(4.2*inch)
    
    rowHeights = (checks_size[0])*[1*inch]

    t=Table(combinedParams,colwidths, rowHeights)
    styled = []
    styled = [ ('TEXTCOLOR',(0,0),(0,-1),colors.blue),
               ('TEXTCOLOR',(0,0),(-1,0),colors.blue),
               ('ALIGN',(0,0),(-1,-1),'CENTER'),
               ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
               ('INNERGRID', (1,1), (-1,-1), 0.25, colors.black),
               ('BOX', (1,1), (-1,-1), 0.25, colors.black),
               ('BACKGROUND',(1,1),(-1,-1),colors.lightgreen)
            ]    

    t.setStyle(TableStyle(styled))
    elements.append(t)

    def coord(x, y, unit=1):
        x, y = x * unit, height -  y * unit
        return x, y

    c = canvas.Canvas("a.pdf", pagesize=letter)
    t.wrapOn(c, width, height)
    t.drawOn(c, *coord(1.8, 9.6, cm))
    c.save()



    # # The order of this setvals is a bit weird, so reverse the set, and turn it into a list:
    # headings = natural_sort(list(setvals))
    # headings.insert(0,'')
    
    # checks =[]
    # checks.append(headings)
    # for p in pats:
    #     curr=[]
    #     curr.append(p)
    
    #     #Recall the first one is blank to match up with the patients column
    #     for k in headings[1:]:
            
    #         if master_list[p][k][0] == '':
    #             curr.append('No')
    #         else:
    #             curr.append('Yes')
    #     checks.append(curr)



    # write the document to disk
    doc.build(elements)             
