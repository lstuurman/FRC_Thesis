// First attempt at script in JS
// Imports :

let CPM = require('../../../cpmjs/build/cpm-cjs.js')
fs = require("fs")
//let FRC = require('../img/frc.json')


// SETUP fn

function setup_sim(l_dir,persist){

    let config =  {
        ndim : 3,
        field_size : [50,50,50],
        conf : {
            torus : [true,true,true],
            T : 20,
            nCellKinds : 2,

            // Adhesion:
            J : [[0,0,5],
                [0,0,-5],
                [5,-5,10]],
            // Volume: 
            V : [0,0,1000],
            LAMBDA_V : [0,0,25],
            // Perimeter:
            P : [0,0,5400],
            LAMDA_P : [0,0,.2],
            // Pref Dir:
            PERSIST : [0,0,persist],
            LAMBDA_DIR : [0,0,l_dir],
            DELTA_T : [0,0,15]

        },
        simsettings : {
	
            // Cells on the grid
            NRCELLS : [0,125],					// Number of cells to seed for all
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
    var fname = '../../data/cpmjs/prefdir/' + String(l_dir) + '_' + String(persist) + 'log.txt'
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

// Costum grid initializer to seed frc:

function initializeGrid(){
		
    // add the initializer if not already there
    if( !this.helpClasses["gm"] ){ this.addGridManipulator(); }

    let nrcells = this.conf["NRCELLS"], cellkind, i;
    let FRC = require('../img/frc.json')
    //console.log(FRC)
    // seed frc cells : 
    for (let x = 0; x < 50; x ++){
        for (let y = 0; y < 50; y++){
            for (let z = 0;z < 50; z++){
                // console.log(x,y,z)
                // console.log(FRC[x][y][z])
                if (FRC[x][y][z] == 1){
                    this.gm.seedCellAt(1,[x,y,z])
                }
            }
        }
    }

    
    // Seed the right number of cells for each cellkind
    for( cellkind = 1; cellkind < nrcells.length; cellkind ++ ){
        
        for( i = 0; i < nrcells[cellkind]; i++ ){
            // first cell always at the midpoint. Any other cells
            // randomly.				
            
            this.gm.seedCell( cellkind+1 );
            
        }
    }
}


// Actual simulations : 
function paramSearch(){
    let persist = Array.from(Array(21).keys(), x => x/20)
    for (let p of persist){
	console.log(p)
        setup_sim(2000,p)
    }
}

paramSearch()
// trial run : 
// setup_sim(2000,.7)


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



    // // Add constraints : 
    // C.add(new CPM.Adhesion(C.conf))
    // C.add(new CPM.VolumeConstraint(C.conf))
    // C.add(new CPM.PerimeterConstraint(C.conf))
    // C.add(new CMP.PersistencConstraint(C.conf))

    // let C_cnv = new CPM.Canvas(C,{zoom : 1})
