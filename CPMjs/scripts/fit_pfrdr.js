let CPM = require('../../../cpmjs/build/cpm-cjs.js')
fs = require("fs")

// SCRIPT FOR FITTING LAMBDA AND PERISTIST OF PREFERENTIAL DIRECTION MODEL

function setup_sim(PERS,LAMBDA){

    let config =  {
        ndim : 3,
        field_size : [64,64,64],
        conf : {
            torus : [true,true,true],
            T : 70,
            nCellKinds : 1,

            // Adhesion:
            J: [[0, 0], [0, 100]],
            // Volume: 
            V : [0,150],
            LAMBDA_V : [0,250],
            // Perimeter:
            P : [0,1400],
            LAMDA_P : [0,2],
            // Pref Dir:
            PERSIST : [0,PERS],
            LAMBDA_DIR : [0,LAMBDA],
            DELTA_T : [0,15]

        },
        simsettings : {
	
            // Cells on the grid
            NRCELLS : [parseInt(64*64*64/150)],					// Number of cells to seed for all
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
    console.log('PERSIST : ',PERS)
    var fname = '../../data/cpmjs/full_PRFDR/Lambda' + String(LAMBDA) + 'PERS' +String(PERS) +'log.txt'
    // create empty file to append to :
    let new_file = fs.writeFile(fname,'',function (err) {
        if (err) throw err;
        console.log('Saved cell tracks');
      })
    logger = fs.createWriteStream(fname, {
        flags: 'a' // 'a' means appending (old data will be preserved)
    })

    let custommethods = {
        // initializeGrid : initializeGrid,
        logStats : logStats
    }

    let sim = new CPM.Simulation( config, custommethods )
    sim.run()
    logger.end()

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
        let line = String(this.time) + "\t" + String(cid) + "\t" + 
        String(this.C.cellKind(cid)) + "\t" + thecentroid.join("\t") + "\n"
        //console.log(line)
        logger.write(line);
        
    }
} 

function fit_PRFDR(){
    let persists = Array.from(Array(11).keys(), x => x/10)
    let Lambdas = Array(500,750,1000,2000,3000,4000,5000)
    // Loop over all combinations : 

    for(var pers of persists){
	if(pers == 0){
	    pers = 0.01
	}
        for(var lamb of Lambdas){
            //console.log(pers,lamb)
            setup_sim(pers,lamb)
        }
    }
}

fit_PRFDR()

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
