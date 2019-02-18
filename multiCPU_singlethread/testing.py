for i in range(0,735):
    for j in range(0,504):
        if abs(data2[j][i] - data1[j][i]) <= 35:
            resized[j][i] = int((data1[j][i]+ data2[j][i] )/2)
        else:
            resized[j][i] = 0

for i in range(0,735):
    for j in range(0,504):
        data3[j][i] = resized[j][i]
        if resized[j][i] == 0:
            k = 0
            m = 0
            W = 1
            for m in range(j-W, j+W):
                for k in range(i-W, i+W):
                    if m > -1 and m < 504 and k > -1 and k < 735:
                        if(resized[m][k] != 0):
                            resized[j][i] = resized[m][k]
                            m  = j+W
                            k = i+W
                            
            if(resized[j][i] ==0):
                print(str(j) + "    " + str(i))
kernel = np.ones((2,2),np.float32)/4
dst = cv2.filter2D(resized,-1,kernel)
cv2.imshow("dst", dst)
fig=plt.figure(figsize=(6, 6))
fig.add_subplot(3, 1, 1)
plt.imshow(data1, cmap = 'gray')
fig.add_subplot(3, 1, 2)
plt.imshow(data2, cmap = 'gray')
fig.add_subplot(3, 1, 3)
plt.imshow(data3, cmap = 'gray')
cv2.imshow("python result", resized)
plt.show()
cv2.waitKey(0)
cv2.destroyAllWindows
