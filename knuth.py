import numpy as np
from numba import cuda
from numba import int16
import time
import sys

# to remove in future
import warnings
from numba.core.errors import NumbaPerformanceWarning
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

BLOCK_SIZE = 32
N = 19456
MAXV = 32767
    
@cuda.jit
def myKernel(c, w, n, c0):
    c1 = cuda.blockIdx.x + c0
    bb = int16(BLOCK_SIZE)
    
    row = cuda.threadIdx.y
    col = cuda.threadIdx.x
    
    C = cuda.shared.array((BLOCK_SIZE, BLOCK_SIZE), dtype=np.int16)
   # A_elements = cuda.shared.array((BLOCK_SIZE, BLOCK_SIZE), dtype=np.int16)
   # B_elements = cuda.shared.array((BLOCK_SIZE, BLOCK_SIZE), dtype=np.int16)
    
    C[row, col] = MAXV
    
    c1_max = min( N // bb, (N + 2 * c0 - 2) // bb)
    if c1 <= c1_max:
        # almost tileable
        for c3 in range(-c0 + c1 + 1, c1):
            
            if row < bb and col < bb:
            
                base4 = bb * (c1 - c0)
                base5 = bb * c1

                c4 = base4 + row
                c5 = base5 + col
                
                if (c4 < 1 or c4 < base4 or c4 > base4 + bb - 1 or c5 < base5 or c5 > base5 + bb - 1 or c5 > N):
                    return
                
                Cvalue = MAXV
                
                cuda.syncthreads()
                
                for c6 in range(bb * c3, bb * c3 + bb):
                    Cvalue = min(Cvalue, w[c4][c5] + c[c4][c6] + c[c6][c5])

                cuda.syncthreads()
                   
                old = C[row, col]
                if Cvalue < old:
                    C[row, col] = Cvalue
                    
        if c0 > 1:
            c[(-c0 + c1) * bb + bb - 1][c1 * bb] = C[bb - 1][0]
                    
        # problematic
        start = max(2, bb * c0 - (bb - 2))
        end   = min(min(min(N - 1, bb * c0 + bb-1), N + bb * c0 - bb * c1), bb * c1 + bb - 2)

        for c4 in range(start, end + 1):
            if c0 >= 1:
                lb = max(max(bb * c1, -bb * c0 + bb * c1 + c4), c4 + 1)
                ub = min(min(N, bb * c1 + bb-1), -bb * c0 + bb * c1 + c4 + bb-1)
                
                c5 = cuda.threadIdx.x + lb
                
                if c5 <= ub:
                    z = int16(c[-c4 + c5, c5])
                    
                    if cuda.threadIdx.y == 0:
                        for c6 in range(-c4 + c5 + 1, -bb * c0 + bb * c1 + bb):
                            z = min(z, w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5])
                            
                        z = min(C[(-c4 + c5) % bb][c5 % bb], c[-c4 + c5][c5])
                        
                        for c6 in range(bb * c1, c5):
                            z = min(z, w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5])
                        
                        c[-c4 + c5, c5] = z                     
            else:
            
                lb = max(bb * c1 + c4, c4 + 1)
                ub = min(N, bb * c1 + bb-1)
                
                c5 = cuda.threadIdx.x + lb
                
                if c5 <= ub:
                    z = int16(c[-c4 + c5, c5])

                    for c6 in range(-c4 + c5 + 1, c5):
                        z = min(z, w[-c4 + c5][c5] + c[-c4 + c5][c6] + c[c6][c5])
                        
                    c[-c4 + c5, c5] = z

# END ------------------------------------------------------------------------------------------------------
                        
print("cuda available:", cuda.is_available())
print("gpus:", list(cuda.gpus))

c = np.zeros((N, N), dtype=np.int16)
w = np.zeros((N, N), dtype=np.int16)
c_CPU = np.zeros((N, N), dtype=np.int16)


t0 = time.perf_counter()
# kopiowanie do GPU
dc = cuda.to_device(c)
dw = cuda.to_device(w)


numBlocks = N // BLOCK_SIZE
bb = BLOCK_SIZE

print("Blocks: " + str(numBlocks) + ", BB: " + str(BLOCK_SIZE))

dimBlock = (BLOCK_SIZE, BLOCK_SIZE)

for c0 in range((N - 1) // bb + 1):
    numBlocks = min((N - 1) // bb, (N + 2*c0 - 2) // bb) - c0 + 1

    myKernel[(numBlocks,), dimBlock](dc, dw, N, c0)
    cuda.synchronize()
    
dc.copy_to_host(c)

t1 = time.perf_counter()
print("czas [s] =", t1 - t0)



'''
for i in range(N - 1, -1, -1):
    for j in range(i + 1, N):
        for k in range(0, j - i):
            S_CPU[i, j] = max(S_CPU[i, k + i] + S_CPU[k + i + 1, j], S_CPU[i, j])

        S_CPU[i, j] = max(S_CPU[i, j], S_CPU[i + 1, j - 1] + paired_host_ascii(seq[i], seq[j]))


for i in range(N):
    for j in range(N):
        if S[i, j] != S_CPU[i, j]:
            print(i, j, ":", S[i, j], S_CPU[i, j])
            print("error")
            raise SystemExit(1)
'''
