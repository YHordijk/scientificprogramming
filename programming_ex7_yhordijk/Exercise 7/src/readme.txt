Implemented algorithm: Shellsort since it is fast for many cases.

We see that for some test-cases involving duplicates that the shellsort algorithm is unstable.
Stability is defined as 

i < j and A[i] = A[j] and m < n

for some array A with with two elements starting in i and j, which are then moves to indices
m and n. This relation does not necessarily hold for shellsort, so it is unstable as is seen in
test-case defined in test004. If this is important a different sorting algorithm should be 
implemented, it could be important for multi-dimensional arrays, such as in the exercise.
It is therefore not very bad for test004 to fail.