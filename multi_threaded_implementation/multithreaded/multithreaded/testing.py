import cv2
import numpy as np
fileName = "C:\\Users\\moham\\Desktop\\multicpu_project\\single_threaded\\images\\im0.png"
img = cv2.imread(fileName,cv2.IMREAD_UNCHANGED)
temp = cv2.resize(img, (735,504))
resized = cv2.cvtColor(temp, cv2.COLOR_BGR2GRAY)

for i in range(0,1):
    for j in range(0,5):
        #resized[j][i] = int(img[j*4][i*4][2]*0.2126 + img[j*4][i*4][1]*0.7152 + img[j*4][i*4][0]*0.0722)
        print(img[i][j])

