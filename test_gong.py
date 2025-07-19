# -*- coding: utf-8 -*-
"""
Created on Fri May 16 17:21:08 2025

@author: valer
"""

import cv2 as cv2
import numpy as np
import matplotlib.pyplot as plt
import Inti_functions as itf
from skimage.metrics import structural_similarity as ssim


def sift_similarity (img1,img2):
    # Initialiser SIFT
    sift = cv2.SIFT_create()
    
    # Détecter les points clés et descripteurs
    kp1, des1 = sift.detectAndCompute(img1, None)
    kp2, des2 = sift.detectAndCompute(img2, None)
    
    # Utiliser FLANN pour matcher les descripteurs
    FLANN_INDEX_KDTREE = 1
    index_params = dict(algorithm=FLANN_INDEX_KDTREE, trees=5)
    search_params = dict(checks=50)
    
    flann = cv2.FlannBasedMatcher(index_params, search_params)
    
    # Trouver les correspondances avec KNN
    matches = flann.knnMatch(des1, des2, k=2)
    
    # Appliquer le ratio test de David Lowe
    good_matches = []
    for m, n in matches:
        if m.distance < 0.7 * n.distance:
            good_matches.append(m)
    
    # Afficher les correspondances
    img_matches = cv2.drawMatches(
        img1, kp1, img2, kp2, good_matches, None,
        flags=cv2.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS
    )
    
    plt.figure(figsize=(12, 6))
    plt.imshow(img_matches)
    plt.title(f"Nombre de bons matches : {len(good_matches)}")
    plt.axis('off')
    plt.show()
    
    # Critère simple : considérer les images similaires si assez de bons matchs
    if len(good_matches) > 10:
        print("✅ Les images sont similaires.")
    else:
        print("❌ Les images ne sont pas similaires.")


def orb_similarity (img1, img2) :
      
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


img1_path= "C:\\Users\\valer\\Desktop\\SharpCap Captures\\test gong\\_11_28_58_gong.jpg"
#img1_path= "C:\\Users\\valer\\Desktop\\SharpCap Captures\\test gong\\2025-03-06_gong.jpg"
img2_path ="C:\\Users\\valer\\Desktop\\SharpCap Captures\\test gong\\_11_28_58_disk_NS.png"
#img2_path ="C:\\Users\\valer\\Desktop\\SharpCap Captures\\test gong\\_2025_03_06-14_28_15-scan_clahe.png"


img1 = cv2.imread(img1_path, cv2.IMREAD_GRAYSCALE)
img1[img1>180] = 180
img1[img1<10]=0
img1= (((img1-10)/(180-10))*255).astype(np.uint8)

img2a = cv2.imread(img2_path, cv2.IMREAD_UNCHANGED) 


img2a[img2a>65000]=65000
img2=(img2a / 256).astype(np.uint8)

plt.imshow(img2, cmap="grey")
plt.show()


axis=0 
offset=0
flag_disk=True
ih1,iw1=img1.shape
m1= ih1//2

ih2,iw2=img2.shape
m2=ih2//2

img1_c =img1[m1-ih1//10:m1+ih1//10,50:-50]
img2_c =img2[m2-ih2//10:m2+ih2//10,5:-5]

a1,a2 = itf.detect_bord (img1_c, axis, offset, flag_disk)
b1,b2 = itf.detect_bord (img2_c, axis, offset, flag_disk)
print(a1,a2)
print(b1,b2)
ratio = (b2-b1)/(a2-a1)
print(ratio)
img_gong=cv2.resize(img1,(int(ih1*ratio),int(iw1*ratio)), cv2.INTER_AREA)
h,w=img_gong.shape
new_h, new_w = img2.shape
top =(new_h-h)//2
bottom = new_h-h-top
left=(new_w-w)//2
right=new_w-w-left
img_gong= cv2.copyMakeBorder(img_gong, top,bottom,left, right, borderType=cv2.BORDER_CONSTANT, value=1)
plt.imshow(img_gong, cmap="grey")
plt.show()

rayon=(b2-b1)//2
cote=int(rayon*(2**0.5))
mi_cote=cote//2
c= new_h//2
# carré inscrit est c-mi_cote
img2=img2[c-mi_cote:c+mi_cote,c-mi_cote:c+mi_cote]
img_gong=np.array(img_gong[c-mi_cote:c+mi_cote,c-mi_cote:c+mi_cote], dtype='uint8')
img_g=np.copy(img_gong)




r_lum=np.mean(img2)/np.mean(img_g)
print("r_lum : ",r_lum)

img22=np.array(img2, dtype='uint8')
img_g=np.array(img_g*r_lum, dtype='uint8')

plt.imshow(img_g, cmap="grey")
plt.show()
plt.imshow(img22, cmap="grey")
plt.show()

"""
diff1=img_g-(img22)
moy1=np.mean(diff1)**2
plt.title(str(moy1))
plt.imshow(diff1)
plt.show()

diff2=img_g-(np.fliplr(img22))
moy2=np.mean(diff2)**2
plt.title(str(moy2))
plt.imshow(diff2)
plt.show()

diff2=img_g-np.flipud(img22)
moy2=np.mean(diff2)**2
plt.title(str(moy2))
plt.imshow(diff2)
plt.show()

diff2=img_g-np.flipud((np.fliplr(img22)))
moy2=np.mean(diff2)**2
plt.title(str(moy2))
plt.imshow(diff2)
plt.show()
"""

inv=['None','EW','NS','NS-EW']

score_PSNR =[cv2.PSNR(img_g, img22),cv2.PSNR(img_g, np.fliplr(img22)),
             cv2.PSNR(img_g, np.flipud(img22)),cv2.PSNR(img_g, np.flipud(np.fliplr(img22)))]
print(cv2.PSNR(img_g, img22))
print(cv2.PSNR(img_g, np.fliplr(img22)))
print(cv2.PSNR(img_g, np.flipud(img22)))
print(cv2.PSNR(img_g, np.flipud(np.fliplr(img22))))
i=np.argmax(score_PSNR)
print(inv[i])


score1, delta = ssim(img_g, img22, data_range=255,full=True)

score2, delta = ssim(img_g, np.fliplr(img22), data_range=255,full=True)

score3, delta = ssim(img_g, np.flipud(img22), data_range=255,full=True)

score4, delta = ssim(img_g, np.flipud(np.fliplr(img22)), data_range=255,full=True)

score_ssim=[score1,score2,score3,score4]
i=np.argmax(score_ssim)
print(inv[i])



"""
score, n_matches = orb_similarity(img_g, img22)
print(f"Score de similarité ORB NS: {score:.2f} ({n_matches} correspondances)")

score, n_matches = orb_similarity(img_g, np.fliplr(img22))
print(f"Score de similarité ORB NS: {score:.2f} ({n_matches} correspondances)")

score, n_matches = orb_similarity(img_g, np.flipud(img22))
print(f"Score de similarité ORB NS: {score:.2f} ({n_matches} correspondances)")

score, n_matches = orb_similarity(img_g, np.flipud(np.fliplr(img22)))
print(f"Score de similarité ORB NS: {score:.2f} ({n_matches} correspondances)")

sift_similarity(img_g, img22)
sift_similarity(img_g, np.flipud(np.fliplr(img22)))
"""


