(function (){
  let str = "x y z o\n";
	for(let i=1;i<1000;i++){
    let jitter = 1+(Math.random()-0.5)/10.;
    let jitter2 = 1+(Math.random()-0.5)/10;
    str += i + " " + (Math.cbrt(i)*jitter) + " " + 1/(Math.pow(Math.cbrt(i),2)*3.)*jitter2 + " " + Math.cbrt(i) + "\n"
    jitter = 1+(Math.random()-0.5)/10.;
    jitter2 = 1+(Math.random()-0.5)/10;
    str += i + " " + (Math.cbrt(i)*jitter) + " " + 1/(Math.pow(Math.cbrt(i),2)*3.)*jitter2 + " " + Math.cbrt(i) + "\n"
  }
console.log(str)
})()