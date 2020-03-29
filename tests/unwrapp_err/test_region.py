import numpy as np
import matplotlib.pyplot as plt
from skimage import measure


regions = np.zeros((100,200,5))

print(regions.shape)


regions[10:30,10:30,0] = 1
regions[10:30,40:70,0] = 2
regions[50:80,80:140,0] = 3
regions[85:100,100:150,0] = 4
regions[15:40,150:180,0] = 5

regions[10:30,10:70,1] = 1
regions[50:80,80:140,1] = 2
regions[85:100,100:150,1] = 3
regions[15:40,150:180,1] = 4

regions[10:30,10:70,2] = 1
regions[50:80,80:140,2] = 2
regions[85:100,100:150,2] = 3
regions[15:40,150:180,2] = 4

regions[10:30,10:60,3] = 1
regions[50:80,90:120,3] = 2
regions[85:100,110:140,3] = 3
regions[15:40,160:180,3] = 4

regions[10:30,10:65,4] = 1
regions[50:80,90:125,4] = 2
regions[85:100,110:157,4] = 3
regions[15:43,160:190,4] = 4


mask = np.ones((100,200))
for i in range(5):
    mask = mask*regions[:,:,i]

mask[mask>0] = 1
labels = measure.label(mask, background=0)

plt.subplot(2,3,1)
plt.imshow(regions[:,:,0])

plt.subplot(2,3,2)
plt.imshow(regions[:,:,1])

plt.subplot(2,3,3)
plt.imshow(regions[:,:,2])

plt.subplot(2,3,4)
plt.imshow(regions[:,:,3])

plt.subplot(2,3,5)
plt.imshow(regions[:,:,4])



plt.subplot(2,3,6)
plt.imshow(labels)
plt.title("common mask")
plt.show()



