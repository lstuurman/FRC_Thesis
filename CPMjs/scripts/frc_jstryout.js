
fs = require("fs")
let FRC = require('../img/frc.json')
var text = fs.readFileSync("../img/frc.json").toString('utf-8');
var textByLine = text.split("\n")

//console.log(text)
var index = textByLine.findIndex(val=>val[0] > 0 || val[1] > 0)
console.log(index)
console.log(FRC[0][0][0])

// // Loop over FRC :
// for(let slice of FRC){
//     for (let row of slice){
//         for( let vox of row) {
//             if (vox == 1){
//                 console.log(vox)
//             }
//         }
//     }
// }