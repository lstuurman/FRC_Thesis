// script for comparing effect of density of peristance, speed and ordered migration behaviour
// Imports :

let CPM = require('../../../cpmjs/build/cpm-cjs.js')
fs = require("fs")
//let FRC = require('../img/frc.json')


// SETUP fn

function setup_sim(ncells){

    let config =  {
        ndim : 3,
        field_size : [50,50,50],
        conf : {
            torus : [true,true,true],
            T : 20,
            nCellKinds : 1,

            // Adhesion:
            J: [[0, 0], [0, 10]],
            // Volume: 
            V : [0,150],
            LAMBDA_V : [0,25],
            // Perimeter:
            P : [0,1400],
            LAMDA_P : [0,0,.2],
            // Pref Dir:
            PERSIST : [0,0.001],
            LAMBDA_DIR : [0,200],
            DELTA_T : [0,15]

        },
        simsettings : {
	
            // Cells on the grid
            NRCELLS : [0,ncells],					// Number of cells to seed for all
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
    let dens = ncells * 1000 / 50*50*50 
    var fname = '../../data/cpmjs/density/PD_' + String(dens) + 'log.txt'
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

// Costum logstats function to save stats:

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

function Densities(){
    let densis = Array.from(Array(21).keys(), x => x/20)
    console.log(densis)
    for (let i = 1; i < 21;i++){
        ncells = parseInt((50*50*50/1000) * densis[i]) 
        console.log(ncells)
        setup_sim(ncells)
    }
}

Densities()