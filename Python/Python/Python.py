
import numpy as np

a = [[1, 2, 3],
     [4, 5, 6]
    ]

u,s,vh = np.linalg.svd(a)

print(u)

print(vh)

print(s)