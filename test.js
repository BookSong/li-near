var li = require('./index.js');

//make a matrix
console.info('Making a matrix...');
var m = new li.Matrix(2,2);
m.print();

//make a vector
console.info('Making a vector...');
var v = new li.Vector(2);
v.print();

console.info('Everything seems to be in order.');
