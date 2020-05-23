//let CPM = require('../../../artistoo/artistoo/build/artistoo-cjs.js')
let CPM = require('../../../cpmjs/build/cpm-cjs.js')
fs = require("fs")

// SCRIPT FOR FITTING LAMBDA AND PERISTIST OF PREFERENTIAL DIRECTION MODEL

function setup_sim(MAX,LAMBDA){

    let config =  {
        ndim : 3,
        field_size : [64,64,64],
        conf : {
            torus : [true,true,true],
            T : 20,
            nCellKinds : 2,

            // Adhesion:
            J: [[0, 0, 0],
                [0, 0, -5],
                [0,-5, 10]],
            // Volume: 
            V : [0,0,150],
            LAMBDA_V : [0,0,25],
            // Perimeter:
            P : [0,0,1500],
            LAMDA_P : [0,0,.2],
            // ACT:
            LAMBDA_ACT : [0,0,LAMBDA],
            MAX_ACT : [0,0,MAX],
            ACT_MEAN : "geometric"

        },
        simsettings : {
	
            // Cells on the grid
            NRCELLS :[0,1], // [parseInt(64*64*64/500)],// Number of cells to seed for all
                                                // non-background cellkinds.
            // Runtime etc
            BURNIN : 200,
            RUNTIME : 5000,
            RUNTIME_BROWSER : "Inf",
            
            // Visualization
            CANVASCOLOR : "EEEEEE",
            CELLCOLOR : ["000000"],
            SHOWBORDERS : [true],				// Should cellborders be displayed?
            BORDERCOL : ["666666"],
            zoom : 2,							// zoom in on canvas with this factor.
            
            // Output images
            SAVEIMG : false,						// Should a png image of the grid be saved
                                                // during the simulation?
            IMGFRAMERATE : 1,					// If so, do this every <IMGFRAMERATE> MCS.
            SAVEPATH : "/home/lau/GIT/FRC_Thesis/CPMjs_output/img",	// ... And save the image in this folder.
            EXPNAME : "DirectedMotionLinear",					// Used for the filename of output images.
            
            // Output stats etc
            STATSOUT : { browser: false, node: true }, // Should stats be computed?
            LOGRATE : 10							// Output stats every <LOGRATE> MCS.
    
        }
    }
    console.log('LAMBDA : ',LAMBDA) 
    console.log('MAX ACT : ',MAX)
    var fname = '../../data/cpmjs/singleFRC_ACT_150/Lambda' + String(LAMBDA) + 'ACT' +String(MAX) +'log.txt'
    // create empty file to append to :
    let new_file = fs.writeFile(fname,'',function (err) {
        if (err) throw err;
        console.log('Saved cell tracks');
      })
    logger = fs.createWriteStream(fname, {
        flags: 'a' // 'a' means appending (old data will be preserved)
    })

    let custommethods = {
        initializeGrid : initializeGrid,
        logStats : logStats
    }
    
    // run simulation 5 times : 
    for (let i=0;i<3;i++){
        iter = String(i)
        let sim = new CPM.Simulation( config, custommethods )
        sim.run()
    }

    //let sim = new CPM.Simulation( config, custommethods )
    //sim.run()
    logger.end()

}

// Costum grid initializer to seed frc:

function initializeGrid(){
		
    // add the initializer if not already there
    if( !this.helpClasses["gm"] ){ this.addGridManipulator(); }

    let nrcells = this.conf["NRCELLS"], cellkind, i;
    let FRC = require('../img/frc.json')
    //console.log(FRC)
    // seed frc cells : 
    for (let x = 0; x < 64; x ++){
        for (let y = 0; y < 64; y++){
            for (let z = 0;z < 64; z++){
                // console.log(x,y,z)
                // console.log(FRC[x][y][z])
                if (FRC[x][y][z] == 1){
                    this.gm.seedCellAt(1,[x,y,z])
                }
            }
        }
    }

function logStats(){
    
    // compute centroids for all cells
    let allcentroids; 
    let torus = false;
    for( let d = 0; d < this.C.grid.ndim; d++ ){
        if( this.C.grid.torus[d] ){
            torus = true;
        }
    }
    if( torus ){
        allcentroids = this.C.getStat( CPM.CentroidsWithTorusCorrection );
    } else {
        allcentroids = this.C.getStat( Centroids );
    } 
    
    for( let cid of this.C.cellIDs() ){
    
        let thecentroid = allcentroids[cid];
        
        // eslint-disable-next-line no-console
        let line = String(this.time) + "\t" + iter + "\t" + //String(cid)
        String(this.C.cellKind(cid)) + "\t" + thecentroid.join("\t") + "\n"
        //console.log(line)
        logger.write(line);
        
    }
} 

function fit_ACT(){
    //let persists = Array.from(Array(11).keys(), x => x/10)
    //persists = Array.from(Array(10).keys(), x => x/10)
    let max_acts = Array(20,50,75,100,200)
    let Lambdas = Array(50,100,500,1000,2000,5000,10000,20000)
    // Loop over all combinations : 

    for(var max of max_acts){
	// if(pers == 0){
	//     pers = 0.9999999999;
	// }
        for(var lamb of Lambdas){
            //console.log(pers,lamb)
            setup_sim(max,lamb)
        }
    }
}

fit_ACT()

// console.log()
// for (let i = 0; i < 11;i++){
// if (i == 0){
//     ncells = 1
//     setup_sim(ncells)
//     console.log(ncells)
// } else {
//     ncells = parseInt((50*50*50/150) * densis[i]) 
//     setup_sim(ncells)
// console.log('N cells : ',ncells)
// }
// }
