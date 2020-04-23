
// fs = require("fs")
// let FRC = require('../img/frc.json')
// var text = fs.readFileSync("../img/frc.json").toString('utf-8');
// var textByLine = text.split("\n")

// //console.log(text)
// var index = textByLine.findIndex(val=>val[0] > 0 || val[1] > 0)
// console.log(index)
// console.log(FRC[0][0][0])

// type = "text/javascript"
// src = "https://github.com/lstuurman/FRC_Thesis/raw/master/CPMjs/img/frc.json"
// console.log(src)

function saveData(blob, fileName)
{
    var a = document.createElement("a");
    document.body.appendChild(a);
    a.style = "display: none";

    var url = window.URL.createObjectURL(blob);
    a.href = url;
    a.download = fileName;
    a.click();
    window.URL.revokeObjectURL(url);
}

var xhr = new XMLHttpRequest();
xhr.open("GET", 'https://github.com/lstuurman/FRC_Thesis/raw/master/CPMjs/img/frc.json');
xhr.responseType = "blob";

xhr.onload = function () {
    saveData(this.response, 'filename');
};
xhr.send();

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