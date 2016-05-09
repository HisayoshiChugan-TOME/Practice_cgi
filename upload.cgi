#!/usr/bin/python
# -*- coding: utf-8 -*-

import cgi
import cgitb
import os, sys
import cv2, matplotlib, math
import numpy as np
import matplotlib.pyplot as plt

def calc_sad(a,b):
    # ベクトル化 & signed化
    a_v = a.flatten().astype(np.float64)
    b_v = b.flatten().astype(np.float64)
    return abs(a_v-b_v).sum()
    
def calc_ssd(a,b):
    # ベクトル化 & signed化
    a_v = a.flatten().astype(np.float64)
    b_v = b.flatten().astype(np.float64)
    return np.dot(a_v - b_v, a_v - b_v ) 
    
def calc_zncc(a,b):
    # ベクトル化 & signed化
    size = float(a.size)
    a_v = a.flatten().astype(np.float64)
    b_v = b.flatten().astype(np.float64)
    # using sum and dot    
    return ( size * np.dot(a_v,b_v) - a_v.sum() *  b_v.sum() ) / \
        ( math.sqrt( size * np.dot(a_v,a_v) -  a_v.sum() *  a_v.sum() ) * \
          math.sqrt( size * np.dot(b_v,b_v) -  b_v.sum() *  b_v.sum() ) )
    # not using sum and dot
    #ave_a = a_v.sum() / float(a.size) 
    #ave_b = b_v.sum() / float(b.size) 
    #print ave_a, ave_b
    #sum0=0.0; sum1=0.0; sum2=0.0;
    #for (i,j) in zip(a_v,b_v):
    #    sum0 += (i - ave_a) * (j - ave_b)
    #    sum1 += (i - ave_a) * (i - ave_a)
    #    sum2 += (j - ave_b) * (j - ave_b)
    #print sum0,sum1,sum2
    #return sum0 / math.sqrt( sum1 * sum2 )

def calc_mi(a,b):
    # ベクトル化 & signed化
    a_v = a.flatten().astype(np.float64)
    b_v = b.flatten().astype(np.float64)
    a_hist = np.zeros(256, dtype=np.float64) # 1d histgram
    b_hist = np.zeros(256, dtype=np.float64) # 1d histgram
    ab_hist = np.zeros( (256,256), dtype=np.float64) # 2d histgram
    # create 1d histgram
    #a_hist =  cv2.calcHist([a],[0],None,[256],[0,256])
    #b_hist =  cv2.calcHist([b],[0],None,[256],[0,256])
    #ab_hist = np.histogram2d(a, b, bins=[256,256]) ??????????
    # ? 2d-hsitgram
    #print "calc 1d hist"
    for i in a_v:
        a_hist[ int(i) ] += 1.0 
    for i in b_v:
        b_hist[ int(i) ] += 1.0   
    # create 2d histgram
    for i,j in zip(a_v,b_v):
        ab_hist[ int(i) ][ int(j) ] += 1.0 
    # normalized hist
    inv = 1.0 / float( a_hist.sum() )
    pa = inv * a_hist #map(lambda n:n*inv, items) 
    pb = inv * b_hist 
    pab= inv * ab_hist 
    # calc mi
    value_mi = 0.0
    for (i,p1) in enumerate(pa):
        for (j,p2) in enumerate(pb):
            if p1*p2 != 0 and  pab[i][j] != 0 : 
                value_mi += pab[i][j] * np.log2( pab[i][j]/(p1*p2) )  
    return value_mi

#######################################################################

cgitb.enable()

#sys.setdefaultencoding("utf-8")
#print 'Content-Type: text/html\n'

form = cgi.FieldStorage()
result1 = 'x'

if form.has_key('file1'):
    item = form['file1']
    if item.file:
        fout = file(os.path.join('/tmp', item.filename), 'wb')
        result1 = os.path.join('/tmp', item.filename)
        while True:
            chunk = item.file.read(1000000)
            if not chunk:
                break
            fout.write(chunk)
        fout.close()
        #result1 = "あがったで1"
       
result2 = 'x'
if form.has_key('file2'):
    item = form['file2']
    if item.file:
        fout = file(os.path.join('/tmp', item.filename), 'wb')
        result2 = os.path.join('/tmp', item.filename)
        while True:
            chunk = item.file.read(1000000)
            if not chunk:
                break
            fout.write(chunk)
        fout.close()
        #result2 = "あがったで2"

im1 = cv2.imread(result1,0)
im2 = cv2.imread(result2,0)

# get shape to resize
im_type = im1.dtype  
re_shape = (0,0,0); re_size = (0,0)
if( im1.size < im2.size ):
    re_shape = im1.shape
    re_size = (re_shape[1], re_shape[0]) # x,yの順 すなわち逆
else:
    re_shape = im2.shape
    re_size = (re_shape[1], re_shape[0]) # x,yの順 すなわち逆

# 小さいのがorg, 大きのがreということにしてく
# 毎回im1,im2 を if else で選択するのはアホらしい(ここで終わり)
im_original = np.ndarray( re_shape, im_type  )
im_resize   = np.ndarray( re_shape, im_type  )
if( im1.size < im2.size ):
    
    im_original = im1.copy()
    im_resize = cv2.resize(im2, re_size )
else:
    im_original = im2.copy()
    im_resize = cv2.resize(im1, re_size )


if ( form["text1"].value  == "xml" ):
    print 'Content-type: text/xml\n'
    print '<?xml version="1.0" encoding="utf-8"?>'
    print '<data>'
    print '<formtype>' + form["text1"].value
    print '<SAD>'+str(calc_sad(im_original, im_resize))+'</SAD>'
    print '<SSD>'+str(calc_ssd(im_original, im_resize))+'</SSD>'
    print '<ZNCC>'+str(calc_zncc(im_original, im_resize))+'</ZNCC>'
    print '<MI>'+str(calc_mi(im_original, im_resize))+'</MI>'
    print '</formtype>'
    print '</data>'
elif ( form["text1"].value == "json" ):
    print 'Content-type: application/json; charset=utf-8\n'
    print '{'
    print '    "' + form["text1"].value + '": {'
    print '        ' +'"SAD": "' + str(calc_sad(im_original, im_resize))  + '",'
    print '        ' +'"SSD": "' + str(calc_ssd(im_original, im_resize))  + '",'
    print '        ' +'"ZNCC": "' + str(calc_zncc(im_original, im_resize)) + '",'
    print '        ' +'"MI": "' + str(calc_mi(im_original, im_resize))   + '"'
    print '    ' +'}'
    print '}'
else:      
    print ('Content-type: text/html; charset=utf-8\n')
    #print 'Content-Type: text/html\n'key in dict 
    print ('<html>')
    print ('<head>')
    print ('<title>Similarities between two images.</title>')
    print ('</head>')
    print ('<body>')
    print (result1 + " vs " + result2)
    print ('<br>')
    print "SAD  = ", calc_sad(im_original, im_resize)
    print ('<br>')
    print "SSD  = ",  calc_ssd(im_original, im_resize)
    print ('<br>')
    print "ZNCC = ",  calc_zncc(im_original, im_resize)
    print ('<br>')
    print "MI   = ",  calc_mi(im_original, im_resize)
    print ('</body>')
    print ('</html>')
