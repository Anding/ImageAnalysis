## Display function algorithm

### Utility functions

1. Clipping function
```
    CLIP(x, s, h) =
    0               if x < s
    1               if x > h
    (x-s)/(h-s)     otherwise
```
2. Midtones transfer function
```
    MID(x, m) =
    0                               if x = 0
    1/2                             if x = m
    1                               if x = 1
    ((m-1) * x) / ((2m-1)*x - m))   otherwise    
```   

### Display algorithm
For each pixel value `x` and using image parameters `m, s, and h`

1. compute `xc = CLIP(x, s, h)`

2. compute `xm = MID(xc, m)`
    
3. plot the pixel using value `xm`

### Computation of image parameters

1. define global parameters
```
    B = 0.25
    C = -2.8
```
    
2. measure image statistics
```
    M = the median pixel value computed from the histogram
    MAD = the median of Abs(x - M) for each pixel value x
    MADN = 1.4826 * MAD
```
    
3. compute the image parameters
```
    a = 0   if M <= 1/2
    a = 1   if M > 1/2

    s = 0   if a = 1
    s = 0   if MADN = 0
    s = min(1, max(0, M + C * MADN))
    
    h = 1   if a = 0
    h = 1   if MADN = 0
    h = min(1, max(0, M - C * MADN))
    
    if a = 0 then compute
        t = M - s
        m = MID(t, B)
        
    else if a = 1 then compute
        t = h - M
        m = MID(B, t)
```
