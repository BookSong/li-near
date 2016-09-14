// author li.qi@outlook.com
// 20150716
var colors = require("colors");

var t = Date.now(); // for recording time elapse
var isVerbose = false;

function setVerbose(verbose){
    isVerbose = verbose;
}



function error(msg){
    console.info(('ERROR: ' + msg).bgRed.white + ('  +' + (Date.now() - t) + 'ms').cyan);
    t = Date.now();    
}

function info(msg){
    console.info(('>>' + msg).green + ('  +' + (Date.now() - t) + 'ms').cyan);
    t = Date.now();    
}

function em(msg){
    console.info(('!!' + msg).yellow + ('  +' + (Date.now() - t) + 'ms').cyan);
    t = Date.now();    
}

function verbose(msg){
    if(isVerbose){
        console.info(('- ' + msg).white + ('  +' + (Date.now() - t) + 'ms').cyan);
        t = Date.now();           
    }
}

function warn(msg){
    console.info(('WARNING: ' + msg).bgYellow.white + ('  +' + (Date.now() - t) + 'ms').cyan);
    t = Date.now();
}

// TODO: process.stdout will consume a lot of memory. Still can't find a function equivalent to 'flush'
function Progress(N, interval, forceOutput){
    this.N = N;
    this.interval = interval || 200;
    this.forceOutput = forceOutput;  // if true, print progress even verbose == false;
    this.t = Date.now();
    this.curT = Date.now();
    this.barLen = 40;
    this.spinnerIndex = 0;
    this.spinner = ['\\', '|', '/', '-'];
}
Progress.prototype.show = function(n){
    if(!this.forceOutput && !isVerbose)
        return;
    if(Date.now() - this.curT < this.interval && n != this.N - 1)
        return;
    this.curT += this.interval;

    process.stdout.cursorTo(0);
    process.stdout.clearLine(1);
    // print spinner
    this.spinnerIndex = (this.spinnerIndex + 1) % this.spinner.length;
    process.stdout.write(this.spinner[this.spinnerIndex] + '  ');
    // print time
    process.stdout.write(('@' + (Date.now() - this.t) + 'ms\t').cyan);
    // print progress
    var percent = (n + 1) / this.N;
    process.stdout.write('['.white);
    for(var i = 0; i < this.barLen; i++){
        process.stdout.write(((i / this.barLen) <= percent ? '■' : '□').white);
    }
    process.stdout.write(']\t'.white);
    process.stdout.write((Math.floor(percent * 100) + '%').yellow);
    if(n == this.N - 1)
        process.stdout.write('\r\n');
    // print mem use
    // TODO

    // print ETA
    // TODO
}

module.exports = {
    error : error,
    info : info,
    em : em,
    warn : warn,
    Progress : Progress,
    setVerbose: setVerbose,
    verbose : verbose
};

