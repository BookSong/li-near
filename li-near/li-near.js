

var ffi = require('ffi');
// var debug = require('debug')('linearalgebra');
var colors = require("colors");
var utils = require('./utils');

var Lapack = null;
var LinearAlgebra = {};


// TODO: add libraries on other platforms
    var LIBRARY = {
        darwin : {},
        freebsd : {},
        linux : {x64 : ['/usr/lib64/liblapack.so.3', '/usr/lib64/libblas.so.3']},
        sunos : {},
        win32 : {arm : [], ia32 : [], x64 : ['lapack_win64_MT', 'blas_win64_MT']}
    };


// Define Fortran Data type, lapack was written with Fortran
    var INT = 4;
    var CHAR = 1;
    var FLOAT = 4;
    var DOUBLE = 8;


// auxiliaries, for type conversion
// types prefixed with "F" are Fortran types
    function toInt(f_int){
        return f_int.readInt32LE(0);
    }
    function toFInt(_int){
        var b = new Buffer(INT);
        b.writeInt32LE(_int,0);
        return b;
    }
    function toChar(f_char){
        return String.fromCharCode(f_char.readUInt8(0));
    }
    function toFChar(_char){
        var b = new Buffer(CHAR);
        b.writeUInt8(_char.charCodeAt(0));
        return b;
    }
    function toFloat(f_float){
        return f_float.readFloatLE(0);
    }
    function toFFloat(_float){
        var b = new Buffer(FLOAT);
        b.writeFloatLE(_float);
        return b;
    }
    function toDouble(f_double){
        return f_double.readDoubleLE(0);
    }
    function toFDouble(_double){
        var b = new Buffer(DOUBLE);
        b.writeDoubleLE(_double);
        return b;
    }


// Link to Lapack
    try {
        Lapack = new ffi.Library(LIBRARY[process.platform][process.arch][0], {
            // Singular Value Decomposition
            "sgesvd_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", ]],
            "dgesvd_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", ]],
            // Singular Value  Decompostion, devide-and-conquer algorithm, faster, but not always reliable. Some occations when I tested this, result is wrong. Use this with caution.
            "sgesdd_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", ]],
            "dgesdd_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", ]],
            // computes the inverse of a matrix using the LU factorization computed by SGETRF.
            "sgetri_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            "dgetri_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            // SGETRS solves a system of linear equations    A * X = B  or  A**T * X = B with a general N-by-N matrix A using the LU factorization computed by SGETRF.
            "sgetrs_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            "dgetrs_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            // SGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
            "sgetrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            "dgetrf_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],

        });
        Blas = new ffi.Library(LIBRARY[process.platform][process.arch][1], {
            // Matrix Multiplication
            "sgemm_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            "dgemm_": ["void", ["pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer", "pointer"]],
            // euclidean norm
            "snrm2_": ["float", ["pointer", "pointer", "pointer"]],
            "dnrm2_": ["double", ["pointer", "pointer", "pointer"]],
            // scale a vector
            "sscal_": ["void", ["pointer", "pointer", "pointer", "pointer"]],
            "dscal_": ["void", ["pointer", "pointer", "pointer", "pointer"]]
        });
    } catch (e) {
        utils.error("Not able to open lapack library.");
        // TODO: Add instructions here.
        // utils.error("Error message reads: ", e);
        throw e;
    }


// Vector Class. Fast Vector Class with Buffer.
    function Vector(m,elementSize,fill){
        // NOTE: Vector only allows float and double elements
        if (JsVector.isJsVector(m)) { //clone
            this.m = m.m;
            this.elementSize = elementSize == DOUBLE ? DOUBLE : FLOAT;
            this.data = new Buffer(this.m * this.elementSize);
            for(var i = 0; i < this.m; i++){
                if(!isNaN(m.getAt[i]))
                    this.setAt(i, m.getAt(i));
                else
                    this.setAt(i, 0);
            }
        } else if (Vector.isVector(m)) { //clone
            this.m = m.m;
            this.elementSize = m.elementSize;
            this.data = new Buffer(m.data);
        } else if (Array.isArray(m)) { //construct from js Array
            this.m = m.length;
            this.elementSize = elementSize == DOUBLE ? DOUBLE : FLOAT;
            this.data = new  Buffer(this.m * this.elementSize);
            for(var i = 0; i < this.m; i++){
                if(!isNaN(m[i]))
                    this.setAt(i, m[i]);
                else
                    this.setAt(i, 0);
            }
        } else if (m >= 0) { //new
            this.m = m;
            this.elementSize = elementSize == DOUBLE ? DOUBLE : FLOAT;
            this.data = new Buffer(this.m * this.elementSize);
            if(Buffer.isBuffer(fill)){
                fill.copy(this.data);
            }else if(fill == 'random'){
                for(var i = 0; i < this.m; i++){
                    this.setAt(i, Math.random());
                }
            }
            else if(!isNaN(fill)){
                for(var i = 0; i < this.m; i++){
                    this.setAt(i, fill);
                }
            }else{// fill 0, default.
                this.data.fill(0);
            }
        } else {
            utils.error('Vector >> invalid params: m = ' + m + ', elementSize = ' + elementSize);
            return undefined;
        }
    }
    Vector.isVector = function(i){
        return i instanceof Vector;
    };
    Vector.prototype.setPrecision = function(elementSize){
        if(elementSize != DOUBLE && elementSize != FLOAT){
            utils.error('Vector.setPrecision >> invalid params. elementSize = ' + elementSize);
            return this;
        }
        if(this.elementSize == elementSize)
            return this;
        // else, really do something
        if(elementSize < this.elementSize)
            utils.warn('WARNING!! losing precision.');

        var r = new Vector(this.m, elementSize);
        for(var i = 0; i < this.m; i++){
            r.setAt(i,this.getAt(i));
        }
        this.data = r.data;
        this.elementSize = elementSize;
        return this;
    };
    Vector.prototype.clone = function(){
        return new Vector(this);
    };
    Vector.prototype.getAt = function(i){
        if(i >= this.m){
            utils.warn('Vector.getAt >> i > m. returning 0');
            return 0;
        }
        var op = this.elementSize == DOUBLE ? 'readDoubleLE' : 'readFloatLE';
        return this.data[op](i * this.elementSize);
    };
    Vector.prototype.setAt = function(i,v){
        var op = this.elementSize == DOUBLE ? 'writeDoubleLE' : 'writeFloatLE';
        this.data[op](v, i * this.elementSize);
        return this; // for chaining
    };
    Vector.prototype.subVector = function(start,len){
        if(start < 0 || len <= 0){
            utils.error('Vector.subVector >> invalid params. start = ' + start + ', len = ' + len);
            return undefined;
        }
        var r = new Vector(len, this.elementSize);
        this.data.copy(r.data, 0, start * this.elementSize, Math.min((start + len), this.m) * this.elementSize);
        return r;
    };
    Vector.prototype.scale = function(scalar, option){
        var useBlas = option && option.useBlas;
        if(useBlas){
            var f_n = toFInt(this.m);
            var f_incx = toFInt(1);
            var f_vector = this.data;
            if(this.elementSize == DOUBLE){
                var f_s = toFDouble(scalar);
                Blas.dscal_(f_n,f_s,f_vector,f_incx);
            }else{
                var f_s = toFFloat(scalar);
                Blas.sscal_(f_n,f_s,f_vector,f_incx);
            }
            return new Vector(this.m, this.elementSize, f_vector);
        }else{
            var r = new Vector(this.m, this.elementSize);
            for(var i = 0; i < this.m; i++){
                r.setAt(i, this.getAt(i) * scalar);
            }
            return r;
        }
    };
    Vector.prototype.plus = function(vector){
        if(!Vector.isVector(vector)){
            utils.error("Vector.plus >> invalid vector: " + vector);
            return undefined;
        }
        var m1 = this.m;
        var m2 = vector.m;
        if(m1 != m2)
            utils.error('Vector.add >> m1 != m2. m1 = ' + m1 + ', m2 = ' + m2);
        r = (m1 > m2 ? this.clone() : vector.clone()); // 以前一个vector的精度为准
        for(var i = 0; i < Math.max(m1, m2); i++){
            r.setAt(i, this.getAt(i) + vector.getAt(i));
        }
        return r;
    };
    Vector.prototype.dot = function(vector){//inner product
        if(!Vector.isVector(vector)){
            utils.error('Vector.dot >> invalid params. vector = ' + vector);
            return undefined;
        }
        var m1 = this.m;
        var m2 = vector.m;
        if(m1 != m2){
            utils.error("JsVector.dot >> invalid params. m1 = " + m1 + ", m2 = " + m2);
            return undefined;
        }
        r = 0;
        for(var i = 0; i < m1; i++){
            r += this.getAt(i) * vector.getAt(i);
        }
        return r;
    };
    Vector.prototype.norm1 = function(){//sum of abs(element)
        var r = 0;
        for(var i = 0; i < this.m; i++){
            r += Math.abs(this.getAt(i));
        }
        return r;
    };
    Vector.prototype.norm2 = function(option){//Euclidean distance
        var useBlas = option && option.useBlas;
        if(useBlas){
            var f_vector = this.data;
            var f_n = toFInt(this.m);
            var f_incx = toFInt(1);
            // NOTE: what does incx mean?
            // the blas docs rather vaguely mentioned the meaning of INCX, if at all. http://www.netlib.org/lapack/explore-html/d7/df1/snrm2_8f.html
            // but it seems that it should always be set to 1.
            // e.g. vector:[1,2,3,4], when incx = 0 norm = 0; when incx = 1 norm = sqrt(1+2^2+3^2+4^2) as we expected. when incx = 2. it jumps over 2 and 4, so norm = sqrt(1+3^2)
            // if set to 0, norm always returns 0;
            return this.elementSize == DOUBLE ?
                Blas.dnrm2_(f_n,f_vector,f_incx)
                : Blas.snrm2_(f_n,f_vector,f_incx);
        }
        else
            return Math.sqrt(this.dot(this));
    };
    Vector.prototype.normalize = function(){//scale to norm = 1
        var norm = this.norm2();
        return this.scale(1/norm);
    };
    Vector.prototype.join = function(vector){
        var m1 = this.m;
        var m2 = vector.m;
        var r = new Vector(m1 + m2, this.elementSize);
        for(var i = 0; i < m1 + m2; i++){
            r.setAt(i, i >= m1 ? vector.getAt(i - m1) : this.getAt(i));
        }
        return r;
    };
    Vector.prototype.toDiagonalMatrix = function(m, n){
        var r = new Matrix(m, n, this.elementSize, 0);
        for(var i = 0; i < this.m; i++){
            r.setAt(i, i, this.getAt(i, i));
        }
        return r;
    };
    //utils
    Vector.prototype.save = function(savePath){
        var fs = require('fs');
        fs.writeFileSync(savePath, toFInt(this.m));
        fs.appendFileSync(savePath, toFInt(this.elementSize));
        fs.appendFileSync(savePath, this.data);
        return this;
    };
    Vector.prototype.print = function(option){
        var digits = option && option.digits;
        digits = isNaN(digits) ? 2 : digits;
        var maxM = option && option.maxM || 50;
        var mp = this.m;
        utils.info("PRINT: Vector >> Dimension: " + this.m);
        if(mp > maxM){
            mp = Math.min(maxM,mp);
            utils.warn("only the first " + mp +" elements are printed");
        }

        process.stdout.write("[  ".grey);
        for(var i =0; i <= mp; i++){
            if(this.m == i)
                continue;
            var ele = i == mp ? '..........'.substring(0, digits + 2) : this.getAt(i);
            switch(i%2){ // alternate column color
                case 0:
                    process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele));
                    process.stdout.write("  ");
                    break;
                //case 1:
                default:
                    process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele).gray);
                    process.stdout.write('  ');
            }
        }
        process.stdout.write("]\r\n".grey);
        return this;
    };
    Vector.read = function(path){
        var fs = require('fs');
        var buffer = fs.readFileSync(path);
        var m = buffer.readInt32LE(0);
        var elementSize = buffer.readInt32LE(4);
        return new Vector(m, elementSize, buffer.slice(8));
    };


// Matrix Class. Fast Matrix Class with Buffer.
    function Matrix(m,n,elementSize,fill){
        // NOTE: JsMatrix allows NAN elements, such as strings.
        if(JsMatrix.isJsMatrix(m)){//clone
            this.m = m.m;
            this.n = m.n;
            this.elementSize = FLOAT; // default precision of conversion
            this.data = new Buffer(this.m * this.n * this.elementSize);
            for(var i = 0; i < this.m; i++){
                for(var j = 0; j < this.n; j++){
                    this.setAt(i,j,m.getAt(i,j));
                }
            }
        }else if(JsVector.isJsVector(m)){//clone from vector to vertical matrix
            this.m = m.m;
            this.n = 1;
            this.elementSize = FLOAT; // default precision of conversion
            this.data = new Buffer(this.m * this.n * this.elementSize);
            for (var i = 0; i <  this.m; i++){
                this.setAt(i, 0, m.getAt(i));
            }
        }else if(Matrix.isMatrix(m)){//clone
            this.m = m.m;
            this.n = m.n;
            this.elementSize = m.elementSize;
            this.data = new Buffer(m.data);
        }else if(Vector.isVector(m)){//clone from vector to vertical matrix
            this.m = m.m;
            this.n = 1;
            this.elementSize = m.elementSize;
            this.data = new Buffer(m.data);
        }else if(Array.isArray(m) && Array.isArray(m[0])){//construct from js 2 dim Array
            this.m = m.length;
            this.n = m[0].length;
            this.elementSize = FLOAT;
            this.data = new Buffer(this.m * this.n * this.elementSize);
            for (var i = 0; i < this.m; i++){
                for(var j = 0; j < this.n; j++){
                    this.setAt(i,j,m[i][j] || 0);
                }
            }
        }else if(Array.isArray(m)){// construct from js 1 dim array to m X 1 matrix
            this.m = m.length;
            this.n = 1;
            this.elementSize = FLOAT;
            this.data = new Buffer(this.m * this.elementSize);
            for(var i = 0; i < this.m; i++){
                this.setAt(i, 0, m[i] || 0);
            }
        }else if(m >= 0 && n >= 0){//new
            this.m = m;
            this.n = n;
            this.elementSize = elementSize == DOUBLE ? DOUBLE : FLOAT;
            this.data = new Buffer(this.m * this.n * this.elementSize);
            if(Buffer.isBuffer(fill)){
                fill.copy(this.data);
            }else if(fill == 'random'){
                for(var i = 0; i < m; i++){
                    for(var j = 0; j < n; j++){
                        this.setAt(i,j,Math.random());
                    }
                }
            }else if(fill == 'identity'){
                this.data.fill(0);
                for(var i = 0; i < Math.min(m,n); i++){
                    this.setAt(i,i,1);
                }
            }else if(!isNaN(fill)){
                for(var i = 0; i < m; i++){
                    for(var j = 0; j < n; j++){
                        this.setAt(i, j, fill);
                    }
                }
            }else{// defaults to 0
                this.data.fill(0);
            }
        }else{
            utils.error('Matrix >> invalid params: m = ' + m + ', n = ' + n + ', elementSize = ' + elementSize);
            return undefined;
        }
    }
    Matrix.isMatrix = function(i){
        return i instanceof Matrix;
    };
    Matrix.prototype.setPrecision = function(elementSize){
        if(elementSize != DOUBLE && elementSize != FLOAT){
            utils.error('Matrix.setPrecision >> invalid params. elementSize = ' + elementSize);
            return this;
        }
        if(this.elementSize == elementSize)
            return this;
        // else, really do something
        if(elementSize < this.elementSize)
            utils.warn('WARNING!! losing precision.');

        var r = new Matrix(this.m, this.n, elementSize);
        for(var i = 0; i < this.m; i++){
            for(var j = 0; j < this.n; j++){
                r.setAt(i,j,this.getAt(i,j));
            }
        }
        this.data = r.data;
        this.elementSize = elementSize;

        return this;
    };
    Matrix.prototype.clone = function(){
        return new Matrix(this);
    };
    Matrix.prototype.getAt = function(i, j){
        if(i >= this.m || j >= this.n){
            utils.warn('Matrix.getAt >> invalid params. i = ' + i + ', j = ' + j);
            return 0;
        }
        var op = this.elementSize == DOUBLE ? 'readDoubleLE' : 'readFloatLE';
        return this.data[op]((this.m * j + i) * this.elementSize);
    };
    Matrix.prototype.setAt = function(i, j, v){
        var op = this.elementSize == DOUBLE ? 'writeDoubleLE' : 'writeFloatLE';
        this.data[op](v, (this.m * j + i) * this.elementSize);
    };
    Matrix.prototype.getVector = function(index, direction){
        direction = direction || 'v'; // 'v' vertical vector, 'h' horizontal vector. Defaults to 'v'
        switch(direction){
            case 'h':
                var r = new Vector(this.n, this.elementSize);
                index = index % this.m;
                for(var i = 0; i < this.n; i++){
                    r.setAt(i,this.getAt(index,i));
                }
                break;
            case 'v':
            default:
                var r = new Vector(this.m, this.elementSize);
                this.data.copy(r.data, 0, index * this.m * this.elementSize);
        }
        return r;
    };
    Matrix.prototype.setVector = function(index, vector, direction){
        vector = vector.clone().setPrecision(this.elementSize);
        if(direction == 'h'){
            index = index % this.m;
            for(var j = 0; j < this.n; j++){
                this.setAt(index, j, vector.getAt(j));
            }
        }else{ // if(direction == 'v')
            vector.data.copy(this.data, index * this.m * this.elementSize, 0, Math.min(vector.m, this.m) * this.elementSize);
        }
        return this;
    };
    Matrix.prototype.subMatrix = function(offsetM, offsetN, m, n){
        if(offsetM < 0 || offsetM >= this.m || offsetN < 0 || offsetN >= this.n || m < 0 || m > (this.m - offsetM) || n < 0 || n > (this.n - offsetN)){
            utils.error('Matrix.subMatrix >> invalid params. offsetM = ' + offsetM + ', offsetN = ' + offsetN + ', m = ' + m + ', n = ' + n);
            return undefined;
        }
        var r = new Matrix(m, n, this.elementSize);
        for(var j = offsetN; j < offsetN + n; j++){
            this.data.copy(r.data, (j - offsetN) * m * this.elementSize, (j * this.m + offsetM) * this.elementSize, (j * this.m + offsetM + m) * this.elementSize);
        }
        return r;
    };
    Matrix.prototype.scale = function(scalar){//scale. wraps matrix.dot
        if(isNaN(scalar)){
            utils.error('Matrix.scale >> invalid params. scalar = ' + scalar);
            return undefined;
        }
        var identity = new Matrix(this.m, this.n, this.elementSize, 'identity');
        return this.dot(identity,{
            alpha:scalar
        });
    };
    Matrix.prototype.normalize = function(direction){//normalize defaults to normalizing columns
        var direction = direction || 'v';
        // for performance reasons, if direction == 'h', transpose it, normalize and transpose back.
        if(direction == 'h')
            return this.transpose().normalize().transpose();
        else
            var r = this.clone();

        for(var i = 0; i < r.n; i++){
            var v = r.getVector(i);
            v = v.normalize();
            r.setVector(i, v, 'v');
        }

        return r;
    };
    Matrix.prototype.dot = function(matrix,option){//multiplication
        var m = this.m;
        var k = this.n;
        if(k != matrix.m){
            utils.error("Matrix.dot >> matrix dimension mismatch. Col count of matrix1 must equal Row count of matrix2");
            return undefined;
        }
        var n = matrix.n;
        var alpha = option && option.alpha;
        if(isNaN(alpha))
            alpha = 1;
        var beta = option && option.beta || 0;
        var c = option && option.c;
        // conform precision
        if(this.elementSize != matrix.elementSize)
            matrix = matrix.clone().setPrecision(this.elementSize);
        if(Matrix.isMatrix(option && option.c))
            f_c = c.clone().setPrecision(this.elementSize).data;
        else
            var f_c = new Buffer(m * n * this.elementSize);// pointer to a m*n matrix data, which is the product of MatrixA and MatrixB.
        // params converted to fortran forms, prefixed with f_
        // params. See http://www.netlib.org/lapack/explore-html/d4/de2/sgemm_8f.html
        var f_m = toFInt(m);
        var f_k = toFInt(k);
        var f_n = toFInt(n);
        var f_transa = toFChar('n');// 'n' no transpose before multiplication
        var f_transb = toFChar('n');
        var f_alpha = this.elementSize == FLOAT ? toFFloat(alpha) : toFDouble(alpha);// scalar to be applied after matrix multiplication;
        var f_beta = this.elementSize == FLOAT ? toFFloat(beta) : toFDouble(beta);// scalar to be applied for matrix C. Set to 0 to indicate that no MatrixC will be added to the product of matrixA and matrixB.
        var f_a = this.data;
        var f_lda = toFInt(Math.max(1,m));
        var f_b = matrix.data;
        var f_ldb = toFInt(Math.max(1,k));
        var f_ldc = toFInt(Math.max(1,m));
        // call blas subroutine
        if(this.elementSize == FLOAT)
            Blas.sgemm_(f_transa,f_transb,f_m,f_n,f_k,f_alpha,f_a,f_lda,f_b,f_ldb,f_beta,f_c,f_ldc);
        else // if(this.elementSize == DOUBLE)
            Blas.dgemm_(f_transa,f_transb,f_m,f_n,f_k,f_alpha,f_a,f_lda,f_b,f_ldb,f_beta,f_c,f_ldc);
        return new Matrix(m, n, this.elementSize, f_c);
    };
    Matrix.prototype.plus = function(matrix){//add. wraps matrix.dot
        var m1 = this.m;
        var n1 = this.n;
        var m2 = matrix.m;
        var n2 = matrix.n;
        if(m1 != m2 || n1 != n2){
            utils.error('Matrix.plus >> dimension mismatch. matrices must have the same column count and the same row count');
            return undefined;
        }
        var identity = new Matrix(this.m, this.m, this.elementSize, 'identity');
        return identity.dot(this,{
            beta:1,
            c:matrix
        });
    };
    Matrix.prototype.inverse = function(){
        //Lapack doc reference: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
        var m = this.m;
        var n = this.n;
        var elementSize = this.elementSize;
        if(m != n){
            utils.error('Matrix.inverse >> to invert, m must equal n');
            return undefined;
        }
        var f_m = toFInt(m);
        var f_n = toFInt(m);
        var f_a = this.data;
        var f_lda = toFInt(Math.max(1,n));
        var f_ipiv = new Buffer(Math.min(m, n) * INT);
        var f_info = toFInt(0);

        //performing LU factorization
        if(elementSize == DOUBLE)
            Lapack.dgetrf_(f_m, f_n, f_a, f_lda, f_ipiv, f_info);
        else
            Lapack.sgetrf_(f_m, f_n, f_a, f_lda, f_ipiv, f_info);

        var info = toInt(f_info);
        if (info !== 0){
            utils.error('Matrix.inverse >> LU factorization failed!');
            return undefined;
        }

        var f_work = new Buffer(elementSize);
        var f_lwork = toFInt(-1);

        //1st getri run to determine optimal f_work size
        if(elementSize == DOUBLE)
            Lapack.dgetri_(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info);
        else
            Lapack.sgetri_(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info);

        info = toInt(f_info);
        if (info !== 0){
            utils.error('Matrix.inverse >> lwork calculation failed!');
            return undefined;
        }

        //else, continue
        lwork = f_work.readFloatLE(0);
        f_work = new Buffer(lwork * elementSize);
        f_lwork = toFInt(lwork);

        //run 2nd time to get actual result
        if(elementSize == DOUBLE)
            Lapack.dgetri_(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info);
        else
            Lapack.sgetri_(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info);

        info = toInt(f_info);
        if (info !== 0){
            utils.error('Matrix.inverse >> getri calculation failed!');
            return undefined;
        }

        return new Matrix(m, n, elementSize, f_a);
    };
    Matrix.prototype.inverseDiagonal = function(){//faster inverse algorithm if matrix is diagonal. condition m == n NOT required.
        var r = new Matrix(this.n, this.m, this.elementSize, 0);
        for(var i = 0; i < Math.min(this.m, this.n); i++){
            var ele = this.getAt(i, i);
            r.setAt(i, i, ele === 0 ? 0 : 1 / ele);
        }
        return r;
    };
    Matrix.prototype.svd = function(option){//singular value decomposition
        // lapack reference: http://www.netlib.org/lapack/explore-html/d8/d49/sgesvd_8f.html
        var useSDD = option && option.useSDD;
        // Faster Singlar Value Decomposition, using devide-and-conquer algorithm
        // !! DO NOT USE !! It seldom converges with larger matrices
        // lapack doc reference http://www.netlib.org/lapack/explore-html/d8/d67/sgesdd_8f.html
        var m = this.m;
        var n = this.n;
        var elementSize = this.elementSize;
        var f_m = toFInt(m);
        f_n = toFInt(n);
        f_a = this.clone().data;
        f_lda = toFInt(Math.max(1, m));
        if(useSDD){
            var f_jobz = toFChar('A');
            var f_iwork = new Buffer(INT * Math.min(m, n) * 8);
        }else{ //if(!useSDD)
            var f_jobu = toFChar('A');
            var f_jobvt = toFChar('A');
        }

        var f_s = new Buffer(Math.min(m, n) * elementSize);
        var f_u = new Buffer(m * m * elementSize);
        var f_vt = new Buffer(n * n * elementSize);
        var f_ldu = toFInt(m);
        var f_ldvt = toFInt(n);

        //Calculating Parameters
        var lwork = -1;
        var f_work = new Buffer(elementSize);
        var f_lwork = toFInt(lwork);
        var f_info = new Buffer(INT);

        //run once to determine optimal f_work size
        if(useSDD && this.elementSize == FLOAT)
            Lapack.sgesdd_(f_jobz, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_iwork, f_info);
        else if(useSDD && this.elementSize == DOUBLE)
            Lapack.dgesdd_(f_jobz, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_iwork, f_info);
        else if(!useSDD && this.elementSize == DOUBLE)
            Lapack.dgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_info);
        else  // if(!useSDD && this.elementSize == FLOAT)
            Lapack.sgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_info);

        if(elementSize == DOUBLE)
            lwork = f_work.readDoubleLE(0);
        else// if elementSize == FLOAT
            lwork = f_work.readFloatLE(0);

        f_work = new Buffer(lwork * elementSize);
        f_lwork = toFInt(lwork);

        //run 2nd time to get actual result
        if(useSDD && this.elementSize == FLOAT)
            Lapack.sgesdd_(f_jobz, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_iwork, f_info);
        else if(useSDD && this.elementSize == DOUBLE)
            Lapack.dgesdd_(f_jobz, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_iwork, f_info);
        else if(!useSDD && this.elementSize == DOUBLE)
            Lapack.dgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_info);
        else  // if(!useSDD && this.elementSize == FLOAT)
            Lapack.sgesvd_(f_jobu, f_jobvt, f_m, f_n, f_a, f_lda, f_s, f_u, f_ldu, f_vt, f_ldvt, f_work, f_lwork, f_info);

        var info = toInt(f_info);
        if(info > 0){
            utils.error('Matrix.svd >> SBDSDC did not converge.');
            return undefined;
        }
        else if(info < 0){
            utils.error('Matrix.svd >> gesxd: the ' + info + 'th argument had an illegal value');
            return undefined;
        }

        return {
            u: new Matrix(m, m, elementSize, f_u),
            s: new Vector(Math.min(m, n), elementSize, f_s),
            vt: new Matrix(n, n, elementSize, f_vt)
        };
    };
    Matrix.prototype.transpose = function(){
        var r = new Matrix(this.n, this.m, this.elementSize);
        for(var i = 0; i < this.m; i++){
            for(var j = 0; j < this.n; j++){
                r.setAt(j,i,this.getAt(i,j));
            }
        }
        return r;
    };
    Matrix.prototype.join = function(matrix,direction){
        direction = direction || 'h';
        var m1 = this.m;
        var n1 = this.n;
        var m2 = matrix.m;
        var n2 = matrix.n;
        if(direction == 'v'){
            if(n1 != n2){
                utils.error('Matrix.join >> dimension mismatch. n1 = ' + n1 + ', n2 =' + n2);
                return undefined;
            }
            var r = new Matrix(m1 + m2, n1, this.elementSize);
            for(var i = 0; i < r.m; i++){
                for(var j = 0; j < r.n; j++){
                    r.setAt(i,j,i < m1 ? this.getAt(i, j) : matrix.getAt(i - m1, j));
                }
            }
        }else{
            if(m1 != m2){
                utils.error('Matrix.join >> dimension mismatch. m1 = ' + m1 + ', n1 = ' + n1);
                return undefined;
            }
            var r = new Matrix(m1, n1 + n2, this.elementSize);
            for(var i = 0; i < r.n; i++){
                if(i < n1)
                    r.setVector(i, this.getVector(i));
                else
                    r.setVector(i, matrix.getVector(i - n1));
            }
        }
        return r;
    };
    // utils
    Matrix.read = function(path){
        var fs = require('fs');
        var buffer = fs.readFileSync(path);
        var m = buffer.readInt32LE(0);
        var n = buffer.readInt32LE(4);
        var elementSize = buffer.readInt32LE(8);
        return new Matrix(m, n, elementSize, buffer.slice(12));
    };
    Matrix.prototype.save = function(savePath){
        var fs = require('fs');
        fs.writeFileSync(savePath, toFInt(this.m));
        fs.appendFileSync(savePath, toFInt(this.n));
        fs.appendFileSync(savePath, toFInt(this.elementSize));
        fs.appendFileSync(savePath, this.data);
        return this;
    };
    Matrix.prototype.print = function(option){
        var digits = option && option.digits;
        digits = isNaN(digits) ? 2 : digits;
        var maxM = option && option.maxM || 20;
        var maxN = option && option.maxN || 20;
        var mp = this.m;
        var np = this.n;
        utils.info("PRINT: Matrix >> Dimension: " + this.m + " X " + this.n);
        if(mp > maxM || np > maxN){
            mp = Math.min(maxM,mp);
            np = Math.min(maxN,np);
            utils.warn("only the upper-left " + mp + "X" + np + " block is printed");
        }

        for(var i =0; i <= mp; i++){
            if(this.m == i)
                continue;
            process.stdout.write("[  ".grey);
            for(var j =0; j <= np; j++){
                if(this.n == j)
                    continue;
                var ele = (i == mp || j == np) ? '..........'.substring(0, digits + 2) : this.getAt(i,j);
                switch(j%2){ // alternate column color
                    case 0:
                        process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele));
                        process.stdout.write("  ");
                        break;
                    //case 1:
                    default:
                        process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele).gray);
                        process.stdout.write('  ');
                }
            }
            process.stdout.write("]\r\n".grey);
        }
        return this;//for chaining
    };



/*use of JsVector and JsMatrix is not recommanded.
But they provide some basic functions in case you can't have a viable  Lapack library, or if you want to print matrices with string elements*/
// TODO: JsVector and JsMatrix have not been thoroughly tested.


// JsVector Class. More intuitive but slower vector class. Use Vector for large vectors and better performance.
    function JsVector(m, fill) {
        // NOTE: JsVector allows NAN elements, such as strings.
        if (JsVector.isJsVector(m) || Vector.isVector(m)) { //clone
            this.m = m.m;
            this.data = [];
            for(var i = 0; i < this.m; i++){
                this.data[i] = m.getAt(i);
            }
        } else if (Array.isArray(m)) { //construct from js Array
            this.m = m.length;
            this.data = [];
            for (var i = 0; i < this.m; i++) {
                this.data[i] = m[i];
            }
        } else if (m >= 0) { //new
            this.m = m;
            this.data = [];
            if (fill == 'random') {
                for (var i = 0; i < m; i++) {
                    this.data[i] = Math.random();
                }
            } else if (fill) {
                for (var i = 0; i < m; i++) {
                    this.data[i] = fill;
                }
            } else // defaults to 0
                for (var i = 0; i < m; i++) {
                this.data[i] = 0;
            }
        } else {
            utils.error('JsMatrix >> invalid params: m = ' + m);
            return undefined;
        }
    }
    JsVector.isJsVector = function(i){
        return i instanceof JsVector;
    }
    JsVector.prototype.print = function(option){
        var digits = option && option.digits;
        digits = isNaN(digits) ? 2 : digits;
        var maxM = option && option.maxM || 50;
        var mp = this.m;
        utils.info("PRINT: Vector >> Dimension: " + this.m);
        if(mp > maxM){
            mp = Math.min(maxM,mp);
            utils.warn("only the first " + mp +" elements are printed");
        }

        process.stdout.write("[  ".grey);
        for(var i =0; i <= mp; i++){
            if(this.m == i)
                continue;
            var ele = i == mp ? '..........'.substring(0, digits + 2) : this.getAt(i);
            switch(i%2){ // alternate column color
                case 0:
                    process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele));
                    process.stdout.write("  ");
                    break;
                //case 1:
                default:
                    process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele).gray);
                    process.stdout.write('  ');
            }
        }
        process.stdout.write("]\r\n".grey);
        return this;
    }
    JsVector.prototype.clone = function(){
        return new JsVector(this);
    }
    JsVector.prototype.getAt = function(i){
        if(i >= this.m){
            utils.warn('JsVector.getAt >> invalid params. i = ' + i);
            return 0;
        }
        return this.data[i];
    }
    JsVector.prototype.setAt = function(i,v){
        if(i < 0 || i > this.m)
            utils.error('JsVector.setAt >> invalid params. i = ' + i);
        else
            this.data[i] = v;
        return this;
    }
    JsVector.prototype.subVector = function(start,len){
        var r = [];
        for(var i = 0; i < len; i++){
            r[i] = this.getAt(start + i) || 0;
        }
        return new JsVector(r);
    }
    JsVector.prototype.plus = function(jsVector){
        if(!JsVector.isJsVector(jsVector)){
            utils.error("JsVector.plus >> invalid jsVector: " + jsVector);
            return undefined;
        }
        var m1 = this.m;
        var m2 = jsVector.m;
        if(m1 != m2){
            utils.error("JsVector.add >> invalid params. m1 = " + m1 + ", m2 = " + m2);
            return undefined;
        }
        r = [];
        for(var i = 0; i < m1; i++){
            r[i] = this.getAt(i) + jsVector.getAt(i);
        }
        return new JsVector(r);
    }
    JsVector.prototype.dot = function(jsVector){//inner product
        var m1 = this.m;
        var m2 = jsVector.m;
        if(m1 != m2){
            utils.error("JsVector.dot >> invalid params. m1 = " + m1 + ", m2 = " + m2);
            return undefined;
        }
        r = 0;
        for(var i = 0; i < m1; i++){
            r += this.getAt(i) * jsVector.getAt(i);
        }
        return r;
    }
    JsVector.prototype.norm1 = function(){//sum of abs(element)
        var r = 0;
        for(var i = 0; i < this.m; i++){
            r += Math.abs(this.getAt(i));
        }
        return r;
    }
    JsVector.prototype.norm2 = function(){//Euclidean distance
        return Math.sqrt(this.dot(this));
    }
    JsVector.prototype.scale = function(scalar){
        var r = [];
        for(var i = 0; i < this.m; i++){
            r[i] = this.getAt(i) * scalar;
        }
        return new JsVector(r);
    }
    JsVector.prototype.normalize = function(){//scale to norm = 1
        var norm = this.norm2();
        return this.scale(1/norm);
    }
    JsVector.prototype.join = function(jsVector){
        var m1 = this.m;
        var m2 = jsVector.m;
        var r = [];
        for(var i = 0; i < m1 + m2; i++){
            r[i] = i >= m1 ? jsVector.getAt(i - m1) : this.getAt(i);
        }
        return new JsVector(r);
    }
    JsVector.prototype.save = function(savePath){
        require('fs').writeFileSync(savePath,JSON.stringify(this));
        return this;
    }
    JsVector.read = function(path){
        var o = JSON.parse(require('fs').readFileSync(path,{encoding:'utf8'}));
        var r = new JsVector(o.m);
        r.data = o.data;
        return r;
    }


// JsMatrix Class. More intuitive but slower matrix class. Use Matrix for large matrices and better performance.
    function JsMatrix(m,n,fill){
        // NOTE: JsMatrix allows NAN elements, such as strings.
        if(JsMatrix.isJsMatrix(m)){//clone
            this.m = m.m;
            this.n = m.n;
            this.data = [];
            for(var i = 0; i < this.m; i++){
                this.data.push([]);
                for(var j = 0; j < this.n; j++){
                    this.data[i][j] = m.data[i][j];
                }
            }
        }else if(JsVector.isJsVector(m)){//clone from vector to vertical matrix
            this.m = m.m;
            this.n = 1;
            this.data = [];
            for (var i = 0; i <  this.m; i++){
                this.data.push([]);
                this.data[i][0] = m.getAt(i,j);
            }
        }else if(Matrix.isMatrix(m)){//clone

        }else if(Vector.isVector(m)){//clone

        }else if(Array.isArray(m) && Array.isArray(m[0])){//construct from js 2 dim Array
            this.m = m.length;
            this.n = m[0].length;
            this.data = [];
            for (var i = 0; i < this.m; i++){
                this.data.push([]);
                for(var j = 0; j < this.n; j++){
                    this.data[i][j] = m[i][j] || 0;
                }
            }
        }else if(m >= 0 && n >= 0){//new
            this.m = m;
            this.n = n;
            this.data = [];
            if(fill == 'random'){
                for(var i = 0; i < m; i++){
                    this.data.push([]);
                    for(var j = 0; j < n; j++){
                        this.data[i][j] = Math.random();
                    }
                }
            }else if(fill == 'identity'){
                for(var i = 0; i < m; i++){
                    this.data.push([]);
                    for(var j = 0; j < n; j++){
                        this.data[i][j] = i == j ? 1 : 0;
                    }
                }
            }else if(fill){
                for(var i = 0; i < m; i++){
                    this.data.push([]);
                    for(var j = 0; j < n; j++){
                        this.data[i][j] = fill;
                    }
                }
            }else{// defaults to 0
                for(var i = 0; i < m; i++){
                    this.data.push([]);
                    for(var j = 0; j < n; j++){
                        this.data[i][j] = 0;
                    }
                }
            }
        }else{
            utils.error('JsMatrix >> invalid params: m = ' + m + ' n = ' + n);
            return undefined;
        }
    }
    JsMatrix.isJsMatrix = function(i){
        return i instanceof JsMatrix;
    }

    JsMatrix.prototype.print = function(option){
        var digits = option && option.digits;
        digits = isNaN(digits) ? 2 : digits;
        var maxM = option && option.maxM || 20;
        var maxN = option && option.maxN || 20;
        var mp = this.m;
        var np = this.n;
        utils.info("PRINT: Matrix >> Dimension: " + this.m + " X " + this.n);
        if(mp > maxM || np > maxN){
            mp = Math.min(maxM,mp);
            np = Math.min(maxN,np);
            utils.warn("only the upper-left " + mp + "X" + np + " block is printed");
        }

        for(var i =0; i <= mp; i++){
            if(this.m == i)
                continue;
            process.stdout.write("[  ".grey);
            for(var j =0; j <= np; j++){
                if(this.n == j)
                    continue;
                var ele = (i == mp || j == np) ? '..........'.substring(0, digits + 2) : this.getAt(i,j);
                switch(j%2){ // alternate column color
                    case 0:
                        process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele));
                        process.stdout.write("  ");
                        break;
                    //case 1:
                    default:
                        process.stdout.write((typeof ele == "number" ? ele.toFixed(digits) : ele).gray);
                        process.stdout.write('  ');
                }
            }
            process.stdout.write("]\r\n".grey);
        }
        return this;//for chaining
    }
    JsMatrix.prototype.clone = function(){
        return new JsMatrix(this);
    }
    JsMatrix.prototype.dot = function(jsMatrix){//multiplication
        var m = this.m;
        var k = this.n;
        if(k != jsMatrix.m){
            utils.error("JsMatrix.dot >> matrix dimension mismatch. Col count of matrix1 must equal Row count of matrix2");
            return undefined;
        }
        var n = jsMatrix.n;
        var r = [];
        for(var i = 0; i < m; i++){
            r.push([]);
            for(var j = 0; j < n; j++){
                r[i][j] = this.getVector(i, 'h').dot(jsMatrix.getVector(j,'v'));
            }
        }
        return new JsMatrix(r);
    }
    JsMatrix.prototype.plus = function(jsMatrix){//add
        var m1 = this.m;
        var n1 = this.n;
        var m2 = jsMatrix.m;
        var n2 = jsMatrix.n;
        if(m1 != m2 || n1 != n2){
            utils.error('JsMatrix.plus >> dimension mismatch. matrices must have the same column count and the same row count');
            return undefined;
        }
        var r = [];
        for(var i = 0; i < m1; i++){
            r.push([]);
            for (var j = 0; j < n1; j++){
                r[i][j] = this.getAt(i,j) + jsMatrix.getAt(i,j);
            }
        }
        return new JsMatrix(r);
    }
    JsMatrix.prototype.getAt = function(i,j){
        if(i >= this.m || j >= this.n){
            utils.warn('JsMatrix.getAt >> invalid params. i = ' + i + ', j = ' + j);
            return 0;
        }
        return this.data[i][j];
    }
    JsMatrix.prototype.setAt = function(i,j,v){
        if(i > this.m || j > this.n || i < 0 || j < 0)
            utils.error('JsMatrix.setAt >> invalid params. i = ' + i + ', j = ' + j);
        else
            this.data[i][j] = v;
        return this; // for chaining
    }
    JsMatrix.prototype.getVector = function(index,direction){
        direction = direction || 'v'; // 'v' vertical vector, 'h' horizontal vector. Defaults to 'v'
        var r = [];
        switch(direction){
            case 'h':
                index = index % this.m;
                for(var i = 0; i < this.n; i++){
                    r[i] = this.getAt(index,i);
                }
                break;
            case 'v':
            default:
                index = index % this.n;
                for(var i = 0; i < this.m; i++){
                    r[i] = this.getAt(i,index);
                }
        }
        return new JsVector(r);
    }
    JsMatrix.prototype.transpose = function(){
        var r = [];
        var m = this.m;
        var n = this.n;
        for(var i = 0; i < n; i++){
            r.push([]);
            for(var j = 0; j < m; j++){
                r[i][j] = this.getAt(j,i);
            }
        }
        return new JsMatrix(r);
    }
    JsMatrix.prototype.save = function(savePath){
        require('fs').writeFileSync(savePath,JSON.stringify(this));
        return this;
    }


//输出
    // 数据类型
    LinearAlgebra.FLOAT = FLOAT;
    LinearAlgebra.DOUBLE = DOUBLE;
    LinearAlgebra.INT = INT;
    LinearAlgebra.CHAR = CHAR;
    // Classes
    LinearAlgebra.JsMatrix = JsMatrix;
    LinearAlgebra.JsVector = JsVector;
    LinearAlgebra.Matrix = Matrix;
    LinearAlgebra.Vector = Vector;

    module.exports = LinearAlgebra;
