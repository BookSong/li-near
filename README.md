# li-near
A easy-to-use linear algebra library with nodejs, using lapack for good performance. Inspired by https://github.com/NaturalNode/node-lapack

# How to use

## install

<code>npm install li-near</code>

### Notes for windows users
If you use x64 system, you can copy contents from ./lapack directory to ~/system32/ or place it in the root directory of your node project. Or you can download your own versions of lapack_win64_MT.dll and blas_win64_MT. Both files are needed. If you use x32 system, you can download x32 versions of the two files and do the same. I haven't figured out a more out-of-the-box way to do this, please feel free to let me know if you have ideas.

### Notes for Linux users
You have to install lapack first. With most Linux distribtutions it is already installed in /usr/lib64/.

### Other systems
I do not have access to other systems, so if you have a solution please let me know. You can make it useful by editing the LIBRARY object in /li-near/li-near.js.

## Getting started
<code>
  var li = require('li-near');
</code>

## Main Objects
### Matrix
Properties:
  m: rows
  n: columns
  elementSize: precision level, can be 4(FLOAT) or 8(DOUBLE)
  data: a buffer object containing the data
  
### Vector
Properties:
  m: rows (Vectors are seen as vertical)
  elementSize: precision level, can be 4(FLOAT) or 8(DOUBLE)
  data: a buffer object containing the data

NOTE: about precision
  When do caculation between matrix and vectors with different precision levels, level will be automatically raised to DOUBLE(high precision).

### JsMatrix/JsVector
Not recommanded to use, but if you want to include other things in a matrix/vector, such as a string, you can try this. I will find time to document their use in the future. But the code is pretty much self-explanatory.

## Make a matrix
### From a 2 dimensional array
<code>
  var matrix = new li.Matrix([[1,2,3],[4,5,6]]);
  matrix.print();
</code>

### From a 1 dimensional array (to a single column matrix)
<code>
  var matrix = new li.Matrix([1,2]);
  matrix.print();
</code>

### From another matrix (clone)
<code>
  var matrix = new li.Matrix(another_matrix);
</code>

### From given parameters: m,n,elementSize,fill
m: row count
n: column count
elementSize: optional, choose from 2 levels of precision. Defaults to FLOAT (single precision), can be set to DOUBLE (double precision), see example below
fill: optional, how to fill the matrix. Can be a buffer object; or 'identity': sets diagonal elements to 1.00 others to 0.00; or 'random': assign random number in the range of [0,1] to each element; or any number, which sets each elements to the same number; if not given, will fill the matrix with 0.00s.
<code>
  var matrix = new li.Matrix(3,4,li.DOUBLE,'identity');
  // li.DOUBLE can be substituted with 8, also works. Any number other than 4 or 8 will be seen as default, which is 4(li.FLOAT).
  matrix.print();
</code>

## Simple calculations
### scale(scalar): scales a matrix by scalar
<code>
  some_matrix.scale(3); // all elements get multiplied by 3.0
</code>

### normalize(direction): normalize vectors in a matrix
direction: 'v'(default) normalize the vertical vectors; 'h' normalize the horizontal vectors;
<code>
  some_matrix.normalize();
</code>

### dot(matrix,option): matrix multiplication
option:
option.alpha, scalar for matrix_A
option.beta, scalar for matrix_B
option.c, another matrix to be added after multiplication, call it matrix_C
The following code does this:
1. scale matrix_A by option.alpha
2. scale matrix_B by option.beta
3. multiply the result from step 1 and 2
4. add matrix_C to it
5. 
option is optional. If omitted, then this does a simple multiplication, as if option is set to:
<code>
{alpha:1,
beta:1,
c:a_zero_matrix}
</code>

<code>
  matrix_A.dot(maxtrix_B,{alpha:2.0, beta:3.0, c: matrix_C);
</code>

### add(matrix): add 2 matrices
<code>
  matrix_A.add(matrix_B);
</code>

### inverse(): get the inverse matrix
<code>
  some_matrix.inverse();
</code>

### inverseDiagonal(): get the inverse matrix, it's faster, but only works on diagonal matrices, whose off-diagonal elements are all zero. It's faster because it only has to calculate the inverse of every diagonal element.
<code>
  some_matrix.inverseDiagonal();
</code>

### transpose(): get transposed matrix
<code>
  some_matrix.transpose();
</code>

## Complex calculations
### svd(): singular value decomposition


## Other functions/ Utilities
### print(): Prints Vector/Matrix/JsVector/JsMatrix.
<code>
  var some_matrix = li.Matrix(3,3);
  var some_vector = li.Vector(3);
  some_matrix.print();
  some_vector.print();
</code>

### read(): reads matrix from file

### save(savePath): save matrix to file

### clone(): clones Vector/Matrix
<code>
  var some_matrix = li.Matrix(3,3);
  var another_matrix = some_matrix.clone();
  // do your magic with another_matrix, some_matrix will not be affected.
</code>

### getAt(m,n): get value at row m and column n
<code>
  some_matrix.getAt(1,1);
</code>

### setAt(m,n,v): set value at row m and column n to v
<code>
  some_matrix.setAt(1,1,0.5);
</code>

### getVector(index, direction): get a Vector from a row/column
index: which row/column do you want. Starting from 0
direction: 'v' (default value) get a column; or 'h' get a row
<code>
  var a_vector = some_matrix(0,'v');
</code>

### setVector(index, vector, direction): set a Vector in a matrix
index: which row/column do you want to change. Starting from 0
vector: the vector you want to put into the matrix
direction: 'v' (default value) set a column; or 'h' set a row

### subMatrix(offsetM, offsetN, m, n): get a subMatrix from a matrix
offsetM: which row to start. Starting from 0
offsetN: which column to start. Starting from 0
m: how many rows to extract
n: how many column to extract

### join(matrix, direction): join 2 matrices

### isMatrix(): returns true if some object is a Matrix
<code>
  li.isMatrix(some_object);
</code>

### setPrecision(elementSize): adjust precision level. elementSize can be 4(li.FLOAT) or 8(li.DOUBLE).
<code>
  some_matrix.setPrecision(li.FLOAT);
</code>

## Supports chaining
<code>
  matrix_D = matrix_A.dot(matrix_B).add(matrix_C).transpose().inverse().print();
  //multiply matrix_A by matrix_B, then add matrix_C to it, then get the transpose of the result, then get the inverse of that, then print the final result, whose value is held in matrix_D.
</code>
