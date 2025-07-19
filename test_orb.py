# -*- coding: utf-8 -*-
"""
Created on Fri May 16 17:21:08 2025

@author: valer
"""

import cv2 as cv2
import numpy as np
import matplotlib.pyplot as plt

def orb_similarity (img1_path, img2_path) :
    
    img1 = cv2.imread(img1_path, cv2.IMREAD_GRAYSCALE)
    img2 = cv2.imread(img2_path, cv2.IMREAD_GRAYSCALE)
    
    orb = cv2.ORB_create(nfeatures=1000)
    kp1, des1 = orb.detectAndCompute(img1,None)
    kp2, des2 = orb.detectAndCompute(img2,None)
    
    bf=cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck= True)
    matches= bf.match(des1, des2)
    matches = sorted(matches, key = lambda x: x.distance)
    
    num_matches=len(matches)
    similarity_score=num_matches/max(len(kp1), len(kp2))
    
    img_matches = cv2.drawMatches(img1, kp1, img2, kp2, matches [:20], None, flags=2)
    plt.figure(figsize=(12,6))
    plt.imshow(img_matches, cmap="grey")
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    
    return similarity_score, num_matches

img1_path= "C:\\Users\\valer\\Desktop\\SharpCap Captures\\2024-02-17\\Capture\\_11_28_58_gong.jpg"
img2_path ="C:\\Users\\valer\\Desktop\\SharpCap Captures\\2024-02-17\\Capture\\_11_28_58_disk_NS.png"
score, n_matches = orb_similarity(img1_path, img2_path)
print(f"Score de similarité ORB NS: {score:.2f} ({n_matches} correspondances)")

img1_path= "C:\\Users\\valer\\Desktop\\SharpCap Captures\\2024-02-17\\Capture\\_11_28_58_gong.jpg"
img2_path ="C:\\Users\\valer\\Desktop\\SharpCap Captures\\2024-02-17\\Capture\\_11_28_58_disk_EW.png"
score, n_matches = orb_similarity(img1_path, img2_path)
print(f"Score de similarité ORB EW: {score:.2f} ({n_matches} correspondances)")


