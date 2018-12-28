/* Copyright 2009 Norbert Stoop, stoopn@ethz.ch */


#include "wire.h"
#include "pthread.h"
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <fstream>


//vtkWireRenderer* GLOBAL_wire_renderer;



void Parameters::setup(std::string filename)
{
std::cout<<"Reading file "<<filename<<"\n";try
  {
    cfg.readFile(filename.c_str());
  }
  catch(const libconfig::FileIOException &fioex)
  {
    std::cout << "I/O error while reading file." << std::endl;
    //return(EXIT_FAILURE);
  }
  catch(const libconfig::ParseException &pex)
  {
    std::cout << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    //return(EXIT_FAILURE);
  }
    std::cout<<"FERTIG\n"<<std::endl;
}


Parameters::Parameters(std::string filename)
{
    setup(filename);

}

void WireSim::read_params()
{
    std::cout<<"Reading parameters\n";
    
	total_sweeps = p.cfg.lookup("SWEEPS");

	torsion_damping_const=p.cfg.lookup("TORSION_DAMPING_CONST");
    position_damping_const=p.cfg.lookup("POSITION_DAMPING_CONST");
    spring_damping_const=p.cfg.lookup("SPRING_DAMPING_CONST");
    bending_damping_const=p.cfg.lookup("BENDING_DAMPING_CONST");
    spring_x0=p.cfg.lookup("SPRING_X0");
    emodul=p.cfg.lookup("EMODUL");
    thickness=p.cfg.lookup("WIRE_DIAMETER");


    max_length=p.cfg.lookup("MAX_LENGTH");//, 0);
    
    dt=p.cfg.lookup("DT");//,0.1);
    fixed_dt=p.cfg.lookup("FIXED_DT");//,false);

    u0x_strain=p.cfg.lookup("U0X_STRAIN");//,0.0);
    u0y_strain=p.cfg.lookup("U0Y_STRAIN");//,0.0);
    u0z_strain=p.cfg.lookup("U0Z_STRAIN");//,0.0);
    u0z_sigma=p.cfg.lookup("U0Z_SIGMA");//,0.0);
    if (p.cfg.exists ("U0Z_INCR_RANDOM"))
    {
        u0z_incr_random=true;
        std::cout<<"Using increment randomness on U0Z\n";
    } else
        u0z_incr_random=false;
    
    if (p.cfg.exists ("EQMODEL"))
    {
        u0x_sigma=p.cfg.lookup("U0X_SIGMA");//,0.0);
        u0y_sigma=p.cfg.lookup("U0Y_SIGMA");//,0.0);
        eqmodel=true;
    } else eqmodel=false;
    
    if (p.cfg.exists ("HAMILTONWALK"))
    {
        hamiltonwalk=true;
        hamiltonu0=1./(spring_x0/sqrt(2.0));
        hamilton_ctr=0;
    } else
        hamiltonwalk=false;

    if (p.cfg.exists ("LOOPMODEL"))
    {
        LM_loopmodel=true;
        LM_avgloopsize=p.cfg.lookup("LM_AVGLOOPSIZE");
        
    } else
        LM_loopmodel=false;
    std::cout<<"Here done\n";

    last_u0z_strain=u0z_strain;
    clamp_nozzle_rotation_r3=p.cfg.lookup("CLAMP_NOZZLE_ROTATION_R3");//,false);
    
    std::cout<<"Here done2\n";

    clamp_free_end_rotation_r3=false; //p.cfg.lookup("CLAMP_FREE_END_ROTATION_R3");//,false);
    
    force_increment_sim=false;
    force_increment_phase=false;
    insertion_force=p.cfg.lookup("INSERTION_FORCE");//,0.0);

    wire_friction_const=p.cfg.lookup("WIRE_FRICTION");//,5.);
    wire_static_friction_const=p.cfg.lookup("WIRE_STATIC_FRICTION");//,15.);
    wall_friction_const=p.cfg.lookup("WALL_FRICTION");//,5.);
    wall_static_friction_const=p.cfg.lookup("WALL_STATIC_FRICTION");//,15.);
    tangential_contact_damping_const=0.;//p.cfg.lookup("TANGENTIAL_CONTACT_DAMPING_CONST");//,0.);

    
    measure_block=p.cfg.lookup("MEASURE_BLOCK");//,10000);
    insertion_speed=p.cfg.lookup("INSERTION_SPEED");//,0.01);
    NMovAv=2*(int)round(spring_x0/(insertion_speed*dt*measure_block));
    initial_deformation=p.cfg.lookup("INITIAL_DEFORMATION");//,0.0);
    xwin_output=p.cfg.lookup("XWIN_OUTPUT");//,true);
    plot_position_block=p.cfg.lookup("PLOT_POSITION_BLOCK");//,1000000);
    qyieldthresh = p.cfg.lookup("QYIELD_THRESHOLD");//,1000);
    prec_last_increase_sweep=0;
    av_epssqrmax.reset();

    
   // debug=parms.value_or_default("DEBUG",false);
    debug=false;
    
    global_friction=p.cfg.lookup("GLOBAL_FRICTION");//, 0.0);
    global_rotational_friction=p.cfg.lookup("GLOBAL_ROTATIONAL_FRICTION");//, global_friction);

    // If Langevin dynamics is used, the global_friction is replaced by the Gamma-parameter of Langevin.
    //
    if (p.cfg.exists ("TEMPERATURE")) {
        langevinsim=true;
        temperature=p.cfg.lookup("TEMPERATURE");
        global_friction=langevin_gamma =p.cfg.lookup("LANGEVIN_GAMMA");
        //langevin_noise = sqrt(2.0*langevin_gamma*temperature/dt); //kb=1 !
        langevin_noise = sqrt(6.0*langevin_gamma*temperature/dt); //kb=1 !
        std::cout<<"WARNING: Performing Langevin dynamics with parameters:\n";
        std::cout<<"T: "<<temperature<<"\nGamma: "<<global_friction<<"\nNoise: "<<langevin_noise<<std::endl;
    } else langevinsim=false;
    
    radius=p.cfg.lookup("RADIUS"); //static_cast<double>(p["RADIUS"]);
    double mass_density=p.cfg.lookup("MASS_DENSITY");//, 1.0);
    
    double border=7*spring_x0;
    block_sweeps=50000;
    radius2=radius*radius;
    container_center.fill(border+radius);
    container_volume = 4./3.*PI*radius2*radius;

    L_x=2.*(radius+border);
    L_y=L_x;
    L_z=L_x;
    total_number_of_nodes=(int)(radius*2/spring_x0)+2;

    sphere_radius=thickness/2.0;
    wire_volume_density = PI*sphere_radius*sphere_radius*spring_x0;

    double area=PI*blitz::pow2(sphere_radius);
    double Itor=0.5*area*blitz::pow2(sphere_radius);
    double Iy= 0.5*Itor + (area*spring_x0)*spring_x0*spring_x0/12.0;
    std::cout<<"Area: "<<area<<", old Iy: "<<Iy;
    Iy=PI*blitz::pow4(sphere_radius)/4.0;
    std::cout<<", new Iy: "<<Iy<<std::endl;

    sphere_mass = area*spring_x0*mass_density;
    sphere_inertia=2.0/5.0*sphere_mass*blitz::pow2(thickness/2.0);


    spring_force_const=area*emodul/spring_x0;
    STICKY_FORCE_CONST=spring_force_const;
    coupling_const=spring_force_const*spring_x0;
    coupling_const*=spring_x0;

    torsion_force_const=emodul*Itor; //We assume the shear modulus G to be equal to the elastic one
    bending_force_const=emodul*Iy;
    if (p.cfg.exists ("BMODUL"))
        bending_force_const=p.cfg.lookup("BMODUL");
    if (p.cfg.exists ("TMODUL"))
        torsion_force_const=p.cfg.lookup("TMODUL");
    if (p.cfg.exists ("FLEXIBLE")) {
        torsion_force_const=0.0;
        bending_force_const=0.0;
        coupling_const=0.0;
        flexible_polymer=true;
    } else
        flexible_polymer=false;
    
    double prefact=1.0;//parms.value_or_default("CONTACT_FORCE_PREFACTOR", 1.0);
    double poisson = 0.33;

    double tEstar= 2.*emodul/(1.-poisson*poisson);
    Estar = prefact*PI/4.*tEstar*spring_x0; //Approximated from two cylinders touching each other.
    wall_force_const=PI/4.*tEstar*spring_x0;

    p.cfg.lookupValue("FILE_PREFIX",file_prefix);
    p.cfg.lookupValue("SEED",seed);
    if (p.cfg.exists ("OUTFILE"))
        p.cfg.lookupValue("OUTFILE",outfile);
    else
        outfile="dump.out";
    if (p.cfg.exists ("INFILE")) {
        p.cfg.lookupValue("INFILE",infile);
        restarted=true;
    }
    //p.cfg.lookupValue("OUTFILE",outfile);
    
    tn2=thickness*2.5;
    tn2sqr=pow(thickness*2.5,2);
    
    LM_u0xcurrent=0.;
    LM_nsegments=0;
    LM_segmentctr=0;
    
}


WireSim::WireSim(std::string filename)
{
std::cout<<"Here\n";
    
sweeps=0;
    vmaxvmax=0.;
    pcorder=4;
    MAXID=1000000;
    block_deform=false;
    all_done=0.;
    rebuild_lc=true;
    exception_occured=false;
    largest_poten=0.;
    plot_block=300; //default pre-init - it's weird but required...
    realsweeptime=0.;
    last_insertion_time=0.;
    restarted=false;
    tnextcol=0;
    old_total_number_of_nodes=total_number_of_nodes;
    during_stiffness_measurement=false;
    stiffness_left=0.;
    stiffness_right=0.;


std::cout<<"There\n";
    p.setup(filename);
    
    read_params();
 	
    srand(seed); 

    elems.resize(0);
    total_number_of_nodes=0; 
    double posdiffy=0.;
    double posdiffz=0.;
    double val;
    do {

	wire::node n;
	wire::edge e;
    //val=container_center[0]+radius-(double)total_number_of_nodes*spring_x0-2.*spring_x0;
    //val=container_center[0]-(double)total_number_of_nodes*spring_x0-2.*spring_x0;
    val=container_center[0]-radius+2.*spring_x0-(double)total_number_of_nodes*spring_x0;
        
	n.pos(0)=val;
	n.pos(1)=container_center[1]+posdiffy;
	n.pos(2)=container_center[2]+posdiffz;
	std::cout<<"Node at position "<<n.pos<<std::endl;

	n.ID=total_number_of_nodes;
	wireiter it = elems.insert(elems.end(),n,e);
	if (total_number_of_nodes>0) {
	  it-1;
	  setup_beam(it);
	} else {
	  it->Q.set(sqrt(2.)/2.,0,-1./sqrt(2.),0.);
	  it->u0 = Tvect(0.); //No initial strain
	  
	}
	++total_number_of_nodes;
    } while (val>container_center[0]-radius-3*spring_x0);
    //  } while (total_number_of_nodes<4);

    for (wireiter it=elems.begin(); it<elems.end(); it++)
      {
	it->Q.display();
      }
    
    old_total_number_of_nodes=total_number_of_nodes;
    if (initial_deformation>0.)
      {
	std::cout<<"Applying initial deformations\n";
	for (int i=0; i<elems.size(); i++)
	  {
	    double rnd=random()-0.5;
	    double diffy=rnd/fabs(rnd)*initial_deformation/50.;
	    rnd = random()-0.5;
	    double diffz=rnd/fabs(rnd)*initial_deformation/50.;
	    posdiffy+=diffy;
	    posdiffz+=diffz;
	    if (fabs(posdiffy)>initial_deformation) {
		posdiffy-=2*diffy;
	    }
	    if (fabs(posdiffz)>initial_deformation) {
		posdiffz-=2*diffz;
	    }
	    elems[i].pos(1)+=posdiffy;
	    elems[i].pos(2)+=posdiffz;
	  }
      }

    setup_matrices();

    channel_radius=thickness/2.; 
	std::cout<<"Simulation created\n";
	std::cout<<"Simulation parameters:\n spring_force_const: "<<spring_force_const
		 <<"\n torsion_force_const: "<<torsion_force_const
		 <<"\n bending_force_const: "<<bending_force_const
		 <<"\n sphere_mass: "<<sphere_mass
		 <<std::endl;

    start();
}


void WireSim::setup_matrices()
{
  IMat.zeros(); IInvMat.zeros(); B1.zeros(); B2.zeros(); B3.zeros(), B01.zeros(); B02.zeros(); B03.zeros(); 
 //   IMat(0,0) = 1./12.*sphere_mass*(3.*thickness*thickness/4. + spring_x0*spring_x0);
//    IMat(1,1) = 1./12.*sphere_mass*(3.*thickness*thickness/4. + spring_x0*spring_x0);
    //Testing!
    //    IMat(2,2) = 0.5*sphere_mass*thickness*thickness/4.;

    IMat(0,0) = sphere_mass*thickness*thickness/4.;
    IMat(1,1) = sphere_mass*thickness*thickness/4.;
    IMat(2,2) = sphere_mass*thickness*thickness/4.;


    IInvMat(0,0) = 1./IMat(0,0);
    IInvMat(1,1) = 1./IMat(1,1);
    IInvMat(2,2) = 1./IMat(2,2);
  
    B1(0,3) =  1.;
    B1(1,2) =  1.;
    B1(2,1) = -1.;
    B1(3,0) = -1.;

    B2(0,2) = -1.;
    B2(1,3) =  1.;
    B2(2,0) =  1.;
    B2(3,1) = -1.;

    B3(0,1) =  1.;
    B3(1,0) = -1.;
    B3(2,3) =  1.;
    B3(3,2) = -1.;


    B01(0,3) =  1.;
    B01(1,2) = -1.;
    B01(2,1) =  1.;
    B01(3,0) = -1.;

    B02(0,2) =  1.;
    B02(1,3) =  1.;
    B02(2,0) = -1.;
    B02(3,1) = -1.;

    B03(0,1) = -1.;
    B03(1,0) =  1.;
    B03(2,3) =  1.;
    B03(3,2) = -1.;

}

void WireSim::setup_beam(wireiter &it) 
{
  wireiter it2=it-1;
 
  it->Q = it2->Q;
 
  it->u0(0)=u0x_strain;
  it->u0(1)=u0y_strain;
    //it->u0(2)=u0z_strain;
  
  if (LM_loopmodel)
  {
    
    if (LM_segmentctr>=LM_nsegments-1) //Loop done, make a new one
    {
        std::cout<<"Determining new loop size"<<std::endl;

        //1. determine new loop size
        // probability of an anchor point, discrete quantity
        double p=spring_x0/LM_avgloopsize;
        double r=rand01();
        LM_nsegments=1;
        while (r>p) {
            r=rand01();
            LM_nsegments++;
        }
        std::cout<<"New loop of size "<<LM_nsegments<<std::endl;
        //2. set u0x accordingly
        double alpha=2*PI/LM_nsegments;
        double rcurv=spring_x0/2./sin(alpha/2.);
        LM_u0xcurrent=1./rcurv;
        it->u0(0)=LM_u0xcurrent;
        //3. determine how long we should wait until the next loop has to be formed.
        LM_segmentctr=0;
    } else {
        it->u0(0) = LM_u0xcurrent;
        LM_segmentctr++;
    }
    //return;
  }
    
  if (hamiltonwalk)
  {
        if (hamilton_ctr<2)
        {
            it->u0(0)=hamiltonu0;
            it->u0(1)=0.;
            it->u0(2)=0.;
            hamilton_ctr++;
        } else if (hamilton_ctr<4) {
            it->u0(0)=0.0;
            it->u0(1)=hamiltonu0;
            it->u0(2)=0.;
            hamilton_ctr++;
        } else if (hamilton_ctr<6) {
            it->u0(0)=hamiltonu0;
            it->u0(1)=0.;
            it->u0(2)=0.;
            hamilton_ctr++;
        }  else if (hamilton_ctr<8) {
            it->u0(0)=0.;
            it->u0(1)=hamiltonu0;
            it->u0(2)=0.;
            hamilton_ctr=0;
        }
        return;
    }
    
  if (eqmodel)
  {
    if (u0x_sigma>0.) it->u0(0)=normalRandom(0.0,u0x_sigma);
    if (u0y_sigma>0.) it->u0(1)=normalRandom(0.0,u0y_sigma);
  }
    
    if (u0z_sigma>0.)
    {
        if (u0z_incr_random)
        {
            it->u0(2)=last_u0z_strain+normalRandom(0.0,u0z_sigma);
            last_u0z_strain=it->u0(2);
        } else
            it->u0(2)=normalRandom(0.0,u0z_sigma);
        
        last_u0z_strain=it->u0(2);
    } else
        it->u0(2)=u0z_strain;
 }


void WireSim::start() {
    double sphere_radius=thickness/2.0;
    
    st_measure_file = file_prefix+"_st_measurements.dat";
    measure_file = file_prefix+"_measurements.dat";
    
    if (restarted)
        load();
    else {
	  measurefile.open(measure_file.c_str());
      measurefile.close();
      
      measurefile.open(st_measure_file.c_str());
      measurefile.close();
      
    }

    std::cout<<"Simulation started\n";
}

void WireSim::halt() {
  std::cout<<"Halting simulation and stopping graphics thread\n";
    save();

  //  pthread_exit(NULL);
}


void WireSim::load() {
    std::cout<<"Loading dump\n";
    std::ifstream ifs(infile.c_str());
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    mili::bistream<> dump(str);

    dump>>sweeps>>total_number_of_nodes;
    dump>>realsweeptime;
	dump>>insertion_force;
    dump>>dt;
	dump>>st;
    dump>>LM_segmentctr;
    dump>>LM_nsegments;
    dump>>LM_u0xcurrent;
    dump>>hamilton_ctr;
    dump>>last_u0z_strain;
    
    std::cout<<"Iterations done so far: "<<sweeps<<" with currently "<<total_number_of_nodes<<" nodes\n";
    
    elems.resize(total_number_of_nodes);

    for (wireiter it=elems.begin(); it<elems.end(); ++it) {
	wire::edge* edge=&it.right_edge();
	
	dump>>*it;
	if (it<elems.end()-1)
		dump>>*edge;
    }
	
}


void WireSim::save() const {
  std::ofstream ofile(outfile.c_str());
    mili::bostream<> dump;
    dump<<sweeps<<total_number_of_nodes;
    dump<<realsweeptime;
    dump<<insertion_force;
    dump<<dt;
    dump<<st;
    dump<<LM_segmentctr;
    dump<<LM_nsegments;
    dump<<LM_u0xcurrent;
    dump<<hamilton_ctr;
    dump<<last_u0z_strain;
    
    std::vector<wire::node> snodes=elems.get_nodes();
    std::vector<wire::edge> sedges=elems.get_edges();

    for (int i=0; i<total_number_of_nodes;++i) {
		dump<<snodes[i];
		if (i<total_number_of_nodes-1)
			dump<<sedges[i];
    }
    
    ofile<<dump.str();
    ofile.close();
}
 /*
bool WireSim::change_parameter(const std::string& name, const alps::StringValue& value)
{
    if(name=="SWEEPS")
	total_sweeps=static_cast<uint32_t>(value);
	else if (name=="DEBUG")
		debug=static_cast<bool>(value);
    else if (name=="MAX_LENGTH")
	max_length=static_cast<int>(value);
    else if (name=="SPHERE_FORCE_CONST")
	sphere_force_const=static_cast<double>(value);
    else if (name=="MEASURE_STIFFNESS_MIN_KIN")
	st.min_kin=static_cast<double>(value);
	else if (name=="MEASURE_STIFFNESS")
	measure_stiffness=static_cast<bool>(value);
    else if (name=="MEASURE_BLOCK")
	measure_block=static_cast<int>(value);
    else if (name=="PLOT_POSITION_BLOCK")
	plot_position_block=static_cast<int>(value);
    else if (name=="INSERTION_SPEED_FACTOR")
	insertion_speed_factor=static_cast<double>(value);
    else if (name=="MAXIMUM_PACKING")
	max_packing=static_cast<double>(value);
    else if (name=="TANGENTIAL_CONTACT_DAMPING_CONST") 
	tangential_contact_damping_const=static_cast<double>(value);
    else if (name=="INSERTION_FORCE") {
	if (insertion_force_start<=0. || insertion_force_end-insertion_force_start<=0.) {
	    insertion_force=static_cast<double>(value);
	    force_increment_sim=false;
	}
    } else if(name=="DT") 
	dt=static_cast<double>(value);
    else if(name=="PREC_MINIMUM_TOL")
	prec_minimum_tolsqr=blitz::pow2(static_cast<double>(value));
    else if(name=="PREC_MAXIMUM_TOL")
	prec_maximum_tolsqr=blitz::pow2(static_cast<double>(value));
    else if(name=="PREC_MINIMUM_DT")
	prec_minimum_dt=static_cast<double>(value);
    else if(name=="PREC_MAXIMUM_DT")
	prec_maximum_dt=static_cast<double>(value);
    else if(name=="FIXED_DT")
		fixed_dt=static_cast<bool>(value);
	
    else if(name=="NODE_MASS") {
	sphere_mass=static_cast<double>(value);	
	sphere_inertia=2.0/5.0*sphere_mass*blitz::pow2(thickness/2.0);
    }
    else if (name=="BLOCK_IMAGE_PLOT") {
	plot_block=static_cast<int>(value);
    } else if (name=="WALL_FRICTION")
	wall_friction_const=static_cast<double>(value);
    else if (name=="INSERTION_SPEED") {
	insertion_speed=static_cast<double>(value); 
    } else if (name=="INSERTION_FORCE") {
	insertion_force=static_cast<double>(value); 
    } else if (name=="TORSION_DAMPING_CONST") {
	torsion_damping_const=static_cast<double>(value);
    } else if (name=="SPRING_DAMPING_CONST") {	
	spring_damping_const=static_cast<double>(value);
    } else if (name=="BENDING_DAMPING_CONST") {	
	bending_damping_const=static_cast<double>(value);
    } else if (name=="T") 
	temperature=static_cast<double>(value);
    else if (name=="MEASURE_CRITICAL_LENGTH") 
	measure_crit_len=static_cast<bool>(value);
    else if (name=="WIRE_FRICTION")
	wire_friction_const=static_cast<double>(value);
    else if (name=="WIRE_STATIC_FRICTION")
	wire_static_friction_const=static_cast<double>(value); 
    else if (name=="WALL_STATIC_FRICTION")
	wall_static_friction_const=static_cast<double>(value); 
    else if (name=="INITIAL_DEFORMATION")
	initial_deformation=static_cast<double>(value);
    else if (name=="VECTOR_IMAGES")
	return true;
    else if (name=="XWIN_OUTPUT") {
	xwin_output=static_cast<bool>(value);
	return true;
    } else if (name=="GLOBAL_FRICTION")
        global_friction=static_cast<double>(value);
    else if (name=="STRAIN_HARDENING_FACTOR")
	strain_hardening_factor=static_cast<double>(value);
    else
	return false; // cannot do it
    return true; // could do it
}
*/

void WireSim::draw_nodes()
{

}

/*  
void WireSim::save_to_file(const boost::filesystem::path& fnpath) const {
    if (!exception_occured) {
	std::cout<<"Storing normal checkpoint\n";
	PCRun::save_to_file(fnpath);
    } else {
	std::cout<<"An exception has occured. No checkpoint is saved\n";
    }

}
*/

void WireSim::reset_all_measurement_variables() {
	BeamBendEnergy=0.0;
	BeamCouplEnergy=0.;
	BeamElongEnergy=0.0;
	BeamTorsEnergy=0.0;
	Pressure=0.0;
	ContactEnergy=0.0;
	CoordinationNumber=0;
	StaticFrictionEnergy =0.;
	KineticEnergy=0.;

}

void WireSim::reset_all_measurement_averages() {
	 av_BeamBendEnergy.reset();
	 av_BeamElongEnergy.reset();
	 av_BeamTorsEnergy.reset();
	 av_BeamCouplEnergy.reset();
	 av_ContactEnergy.reset();
	 av_InsertionForce.reset();
	 av_Pressure.reset();
	 av_CoordinationNumber.reset();
	 av_StaticFrictionEnergy.reset();
	 av_KineticEnergy.reset();
}


void WireSim::dostep() {
    if (sweeps<total_sweeps && all_done<1.) {

	if ((sweeps%plot_block==0 ) || total_number_of_nodes>=old_total_number_of_nodes+plot_block) { //Draw an image every plot_block nodes
	    old_total_number_of_nodes=total_number_of_nodes;
	    draw_nodes();	    
	}
	reset_all_measurement_variables();
	predictor_corrector_step();
	packing_ratio = total_number_of_nodes*wire_volume_density/container_volume;


	if (sweeps%plot_position_block==0 ) {
        save();
std::string pfilename, cfilename;
          char buffer[50];
          sprintf(buffer, "_pos_%09d",  sweeps);
std::cout<<"test\n";
	  pfilename=file_prefix+std::string(buffer)+".vtp";
          sprintf(buffer, "_ct_%09d",  sweeps);
	  cfilename=file_prefix+std::string(buffer)+".dat";
	  output_vtk(pfilename);
	  ctfile.open(cfilename.c_str());
	  for (wireiter it=elems.begin(); it<elems.end(); ++it) {
	    for (int k=0; k<it->contacts.size(); ++k) {
	      if (it->contacts[k].distsqr<=blitz::pow2(1.1*thickness) ) {
		ctfile<<std::distance(elems.begin(),it)<<"\t"<<it->contacts[k].idx<<"\t"<<it->contacts[k].s1t<<"\t"<<it->contacts[k].s2t
		      <<"\t"<<sqrt(it->contacts[k].distsqr)<<"\n";
	      }
	    }
	  }
	  ctfile.close();
	}

/*
	if (measure_stiffness && (st.during_stiffness_measurement || st.measure_next(packing_ratio)))
	  {
	    st.during_stiffness_measurement = true;
	    insertion_speed = 0.0;
	    
	    if (av_KineticEnergy.avg()<st.min_kin && !st.measuring) {
	      st.measuring = true;
	      st.msweeps=0;
	      std::cout<<"At sweep "<<sweeps<<": starting with stiffness measurement (system relaxed with Ekin="<<av_KineticEnergy.avg()<<"<"<<st.min_kin<<")\n";
	      reset_all_measurement_averages();
	    }
	    if (st.measuring) {
	      st.msweeps++;
	      
	      if (!(st.msweeps%measure_block))
		{
		  //Measurements done
		  do_measurements(st_measure_file,true);
		  
		  st.measuring = false;
		  st.during_stiffness_measurement=false;
		  insertion_speed = parms.value_or_default("INSERTION_SPEED",0.01);
		  st.advance();
		  std::cout<<"Next stiffness measurement at packing density "<<st.next_x()<<"\n";
		} else 
		do_measurements(st_measure_file,false);
	    } else 
	      do_measurements(measure_file, !(sweeps%measure_block));
	    
	  } else */
	  do_measurements(measure_file, !(sweeps%measure_block));
 

    if (sweeps%50==1 || langevinsim)
	    rebuild_lc=true;
	sweeps++;
	realsweeptime+=dt;
	if (max_length>0 && total_number_of_nodes>=max_length) 
	    all_done=1.0;
    } else if (all_done!=1.)
	all_done=1.;  
}


void WireSim::do_measurements(std::string filename,bool write_measurements) {
	wireiter it=elems.end()-1;
	wireiter it2=it-1;
	Tvect force(0.);
	 
	//Loop over 4 nodes in channel to average tensile force:
	for (int i=0; i<3; i++)
	{
		Tvect rij = it->pos - it2->pos;
		Tvect vij = it->vel - it2->vel;
		double len = it2.right_edge().length;
		Tvect uij = rij/len;
		
		double elo=(len-spring_x0);
		force += uij*(spring_force_const*elo+spring_damping_const*blitz::dot(vij,uij));
		it--;
		it2--;
	}	
	force/=3.;

  	
	av_BeamBendEnergy.add(BeamBendEnergy);
	av_BeamElongEnergy.add(BeamElongEnergy);
	av_BeamTorsEnergy.add(BeamTorsEnergy);
	av_BeamCouplEnergy.add(BeamCouplEnergy);
	av_ContactEnergy.add(ContactEnergy);
	av_InsertionForce.add(force[0]);
	av_Pressure.add(Pressure);
	av_CoordinationNumber.add(CoordinationNumber);
	av_StaticFrictionEnergy.add(StaticFrictionEnergy);
	av_KineticEnergy.add(KineticEnergy);
	
	//	if (!(sweeps%measure_block))
	if (write_measurements)
	{ //Write out the measurements
	

	  measurefile.open(filename.c_str(), std::ios::app|std::ios::out);
		if (!measurefile.is_open())
		{
		  std::cout<<"At sweep "<<sweeps<<": Measurement file "<<filename<<" could not be opened!\n";
		  exception_occured=true;
		  all_done=1.;
		  return;	
		}
		double real_packing_ratio=0;
		double wdensity = PI*sphere_radius*sphere_radius;
		Tvect dvec;
		double dvecl;
		for (wireiter it=elems.begin(); it<elems.end()-1; it++)
		  {
		    dvec=it->pos-container_center;
		    dvecl=blitz::dot(dvec,dvec);	
		    if (dvecl<=radius2) 
		      real_packing_ratio += it.right_edge().length; 
		  }
		real_packing_ratio*=wdensity/container_volume;

		measurefile<<realsweeptime<<"\t";
		measurefile<<total_number_of_nodes<<"\t";
		measurefile
            <<av_KineticEnergy.avg()<<"\t"
			<<av_BeamBendEnergy.avg()<<"\t"
			<<av_Pressure.avg()/(4.0*PI*radius*radius)<<"\t"
			<<av_CoordinationNumber.avg()<<"\t"
			<<av_InsertionForce.avg()<<"\t"
			<<av_BeamElongEnergy.avg()<<"\t"
			<<av_BeamTorsEnergy.avg()<<"\t"
			<<av_ContactEnergy.avg()<<"\t"
		    <<av_StaticFrictionEnergy.avg()<<"\t"
			<<packing_ratio<<"\t"
			<<real_packing_ratio<<"\t"
			<<av_BeamCouplEnergy.avg()<<"\t"
		        <<av_InsertionForce.max()<<"\t"
            <<KineticEnergy<<"\t"
            <<KineticEnergy/total_number_of_nodes*2./3.
			<<"\n";
		
		measurefile.close();  
		reset_all_measurement_averages();
	}
}



double WireSim::work_done() const {
    return ((all_done==1.)?(1.):(sweeps/(double)total_sweeps));
}


void WireSim::computeBeamForces()
{
	int i=0;
	double len;
	double elo;
	Tvect rij,vij,uij,spring_force,dvreldrix,dvreldriy,dvreldriz,df;
	Quaternion qq1, Qsum, Qdiff;
	Quaternion qq2;
	double ee,f;
	double rx,ry,rz, r2x,r2y,r2z;
	Tvect firstterm;
	double invspring_x0 = 1./spring_x0;
	double prefact2 = position_damping_const/(spring_x0*spring_x0*spring_x0*spring_x0*spring_x0);
	double dprefact = 4./spring_x0*bending_damping_const;
	double dprefact2 = 4./spring_x0*torsion_damping_const;
	double prefact3;
	

	for (wireiter it=elems.begin(); it<elems.end()-1; it++, i++)
   //     for (wireiter it=elems.end()-4; it<elems.end()-1; it++, i++)

	{
		wireiter it2=it+1;
		wireiter it0=it-1;
	
		
		//Spring force:
		rij = it2->pos - it->pos;
		vij = it2->vel - it->vel;
	        len = it.right_edge().length;
		uij = rij/len;
		
	        elo=(len-spring_x0);
		// This is the force due to elongation of the spring connecting the two nodes
		spring_force = uij*(spring_force_const*elo+spring_damping_const*vij.dot(uij));
	
		it->E_tens += 0.25*spring_force_const*elo*elo;  //Linear exptrapolation (mean) from the edges to the nodes
		it2->E_tens += 0.25*spring_force_const*elo*elo;
			
		it->force+=spring_force;
		it2->force-=spring_force;
		
		BeamElongEnergy += 0.5*spring_force_const*elo*elo;
		
        if (flexible_polymer)
            continue;
        
		
		//Relative velocity damping (all degrees of freedom as long as relative to each other)
		double dd = vij.dot(rij);
		Tvect vrelsc = prefact2*rij*dd;
		
		dvreldrix     = -rij*vij(0);
		dvreldrix(0) -= dd;
		dvreldriy     = -rij*vij(1);
		dvreldriy(1) -= dd;
		dvreldriz     = -rij*vij(2);
		dvreldriz(2) -= dd;

		df(0) = vrelsc.dot(dvreldrix);
		df(1) = vrelsc.dot(dvreldriy);
		df(2) = vrelsc.dot(dvreldriz);

		it->force  -= df;
		it2->force += df;

		//Bending and torsion:
		//Here, we take derivatives of the discrete bending energies by the quaternion state variables
		//Quaternions for edge i are stored in nodes i, so we can traverse until total_number_of_nodes-1
		if (it<elems.end()-2) //Bending compares two quaternions, so we can only do it for the second last element
		  {
		    Qsum  = it->Q+it2->Q;
		    Qdiff = it2->Q-it->Q;
		    //Two bending parts:
		    f = (B1mult(Qsum)).dot(Qdiff)*invspring_x0-it->u0[0];
		    ee = f*f*bending_force_const*spring_x0*0.5;
		    BeamBendEnergy += ee;
		    it2->E_bend += 0.5*ee;  //Linear extrapolation (mean) from the edge to the nodes
		    it->E_bend += 0.5*ee;

		    qq1 = B1mult(Qsum);
		    qq2 = B1mult(Qdiff, true);	
		    it->Qa  = it->Qa  - (qq2-qq1)*(bending_force_const*f);
		    it2->Qa = it2->Qa - (qq2+qq1)*(bending_force_const*f);
		    
		    f = (B2mult(Qsum)).dot(Qdiff)*invspring_x0-it->u0[1];
		    ee =f*f*bending_force_const*spring_x0*0.5;
		    BeamBendEnergy += ee;
		    it2->E_bend += 0.5*ee;
		    it->E_bend += 0.5*ee;

		    qq1 = B2mult(Qsum);
		    qq2 = B2mult(Qdiff, true);
		    it->Qa  = it->Qa  - (qq2-qq1)*(bending_force_const*f);
		    it2->Qa = it2->Qa - (qq2+qq1)*(bending_force_const*f);
		    
		    //One torsion/twist part:
		    f = (B3mult(Qsum)).dot(Qdiff)*invspring_x0-it->u0[2];
		    ee = f*f*torsion_force_const*spring_x0*0.5;
		    BeamTorsEnergy += ee;
		    it2->E_tors += 0.5*ee;
		    it->E_tors += 0.5*ee;

		    qq1 = B3mult(Qsum);
		    qq2 = B3mult(Qdiff, true);
		    if (!clamp_free_end_rotation_r3 || it!=elems.begin())
		      it->Qa  = it->Qa  - (qq2-qq1)*(torsion_force_const*f);
		    it2->Qa = it2->Qa - (qq2+qq1)*(torsion_force_const*f);
		    
		    
		    //Relative rotational damping:
		    f = dprefact*(B01mult(it2->Q).dot(it2->Qv) - B01mult(it->Q).dot(it->Qv));
		    it->Qa  = it->Qa  + B01mult(it->Qv, true)*f;
		    it2->Qa = it2->Qa - B01mult(it2->Qv, true)*f;
		 	
		    f = dprefact*(B02mult(it2->Q).dot(it2->Qv) - B02mult(it->Q).dot(it->Qv));
		    it->Qa  = it->Qa  + B02mult(it->Qv, true)*f;
		    it2->Qa = it2->Qa - B02mult(it2->Qv, true)*f;
		    
		    f = dprefact2*(B03mult(it2->Q).dot(it2->Qv) - B03mult(it->Q).dot(it->Qv));
		    if (!clamp_free_end_rotation_r3 || it!=elems.begin())
		      it->Qa  = it->Qa  + B03mult(it->Qv, true)*f;
		    it2->Qa = it2->Qa - B03mult(it2->Qv, true)*f;	    
		  }

		//The quaternions are coupled to the edges so that d3 is parallel to the edge. This can be done for all quaternions, except the last one (which does not yet have an edge):
		
		rx=it->pos(0);
		ry=it->pos(1);
		rz=it->pos(2);
		r2x=it2->pos(0);
		r2y=it2->pos(1);
		r2z=it2->pos(2);
		firstterm = coupling_const*(uij - it->Q.getD3());  //coupling_const = konst*spring_x0 !
		//	it->Q.display();
		//	std::cout<<"D3: "<<it->Q.getD3()<<std::endl;
		//	std::cout<<"Firstterm: "<<firstterm<<std::endl;

		BeamCouplEnergy += firstterm.dot(firstterm)*0.5;

		prefact3 = 1./(len*len*len);
		Tvect ddrx,ddry,ddrz;
		ddrx(0) = (-r2y*r2y - r2z*r2z + 2.*r2y*ry - ry*ry + 2.*r2z*rz - rz*rz);
		ddrx(1) = rij(0)*rij(1) ;
		ddrx(2) = rij(0)*rij(2);

		ddry(0) = rij(0)*rij(1);
		ddry(1) = (-r2x*r2x - r2z*r2z + 2.*r2x*rx - rx*rx + 2.*r2z*rz - rz*rz);
		ddry(2) = rij(1)*rij(2);	

		ddrz(0) = rij(0)*rij(2);
		ddrz(1) = rij(1)*rij(2);
		ddrz(2) = (-r2x*r2x - r2y*r2y + 2.*r2x*rx - rx*rx + 2.*r2y*ry - ry*ry);
	
		Tvect couplingforce;
		couplingforce(0) = prefact3*firstterm.dot(ddrx);
		couplingforce(1) = prefact3*firstterm.dot(ddry);
		couplingforce(2) = prefact3*firstterm.dot(ddrz);
		it->force -= couplingforce;
		it2->force += couplingforce;

		//Now the coupling torques on the quaternions:
		Vector4 couplingtorque;
		couplingtorque(0) = firstterm.dot(it->Q.getdD3dq0());
		couplingtorque(1) = firstterm.dot(it->Q.getdD3dq1());
		couplingtorque(2) = firstterm.dot(it->Q.getdD3dq2());
		couplingtorque(3) = firstterm.dot(it->Q.getdD3dq3());
		it->Qa = it->Qa + Quaternion(couplingtorque[0],couplingtorque[1],couplingtorque[2],couplingtorque[3]);
	}
}


#define SMALL_NUM  0.000000001

double WireSim::segment_segment_distsqr(Tvect S1P0, Tvect S1P1, Tvect S2P0, Tvect S2P1, Tvect &s1p, Tvect &s2p, double &s1t, double &s2t)
{
    Tvect   u = S1P1 - S1P0;
    Tvect   v = S2P1 - S2P0;
    Tvect   w = S1P0 - S2P0;
    double  a = blitz::dot(u,u);        // always >= 0
    double  b = blitz::dot(u,v);
    double  c = blitz::dot(v,v);        // always >= 0
    double  d = blitz::dot(u,w);
    double  e = blitz::dot(v,w);
    double  D = a*c - b*b;       // always >= 0
    double  sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
    double  tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    Tvect dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
    s1t=sc;
    s2t=tc;
    s1p=S1P0+(sc*u);
    s2p=S2P0+(tc*v);
  
    return blitz::dot(dP,dP);   // return the closest distance squared
}


void WireSim::update_collision_forces() {
     
    int i;
    int j;
    int cidx[3];
    int cidx2[3];
    Tvect v;
    Tvect pn;
    //  Tvect dvec2;
    int sidx;
    int sidx2;   
    double t;
    double distsqr;

    Tvect s1p;  //Contact points on segments S1 nd S2
    Tvect s2p;  //Vector from node to contact point right
    double s1t; // Parametric value of contact point on line segment S1 and S2
    double s2t;
    bool just_rebuilt=false;
  //  bool n2n;
    // Build linked cell structure and find contacts: 
    if (rebuild_lc) {
	//	std::cout<<"Sweep "<<sweeps<<": rebuilding LC\n";
      rebuild_lc=false;
       // just_rebuilt=true;
      //	rc=dt*sqrt(vmaxvmax)*2;
      rc=((rc<spring_x0*2.)?(spring_x0*2.):(rc));
//Save bet:
      rc = spring_x0*3.;
		
      rcrc=rc*rc;
      L_cx=L_x/rc;
      L_cy=L_y/rc;
      L_cz=L_z/rc;
        
        if (L_cx>25) {
            L_cx=L_cy=L_cz=25;
            rc=((double)L_x)/((double)L_cx);
            
        }
      //   std:;cout<<L_cz<<std::endl;
      head.clear();
      linkedcell.clear();
      head.resize(L_cx*L_cy*L_cz,-1);
      linkedcell.resize(elems.size(),-1);	 
	for (int ct=0; ct<elems.size(); ++ct) {
	    cidx[0]=(int)((elems[ct].pos/rc)(0));
	    cidx[1]=(int)((elems[ct].pos/rc)(1));
	    cidx[2]=(int)((elems[ct].pos/rc)(2));		
		
	    elems[ct].in_ct_list=false;
	    if (elems[ct].contacts.size()>0) {
		for (std::vector<wire::contact>::iterator kk=elems[ct].contacts.begin(); kk<elems[ct].contacts.end(); ++kk) {
		    if (!kk->valid) {
				kk=elems[ct].contacts.erase(kk);
				--kk;
		    }	
		}
		elems[ct].ct_idx_map.clear();
		elems[ct].ct_idx_map.resize(0);
		for (int kk=0; kk<elems[ct].contacts.size(); ++kk) {
		    elems[ct].ct_idx_map[elems[ct].contacts[kk].idx]=1;
		}
	    }
	    sidx=cidx[0]*L_cz*L_cy+cidx[1]*L_cz+cidx[2];
	    linkedcell[ct]=head[sidx];
	    head[sidx]=ct;
	}
	contact_list.resize(0);
	
	double ncut=thickness/spring_x0+0.1*thickness;
	wire::contact newcontact;
	
	wireiter it=elems.begin();
	for (cidx[0]=0; cidx[0]<L_cx; cidx[0]++) 
	    for (cidx[1]=0; cidx[1]<L_cy; cidx[1]++) {	 
		for (cidx[2]=0; cidx[2]<L_cz; cidx[2]++) {   
		    sidx=cidx[0]*L_cy*L_cz+cidx[1]*L_cz+cidx[2];
		    i = head[sidx];	    
			while (i>-1) {
		    	for (cidx2[0]=cidx[0]-1; cidx2[0]<=cidx[0]+1; cidx2[0]++) 	
				for (cidx2[1]=cidx[1]-1; cidx2[1]<=cidx[1]+1; cidx2[1]++) {
			    for (cidx2[2]=cidx[2]-1; cidx2[2]<=cidx[2]+1; cidx2[2]++) {
					sidx2=((cidx2[0]+L_cx)%L_cx)*L_cz*L_cy+((cidx2[1]+L_cy)%L_cy)*L_cz+((cidx2[2]+L_cz)%L_cz);		    
			
				    j=head[sidx2];
				    
				    while (j>-1) { //We can check symmetrically, since segment-segment is a symmetric interaction
				      // We do not need to check for the directly connected neighbours			
				      if (j>i && i<elems.size()-1 && j<elems.size()-1 && i+1<j && (i<floor(j-ncut) || i>ceil(j+ncut+1))) {
					    v=elems[i].pos-elems[j].pos;
					    					    
					    // This rules out those neighbours that are too far away for sure
					    if ( blitz::dot(v,v)<blitz::pow2(tn2+spring_x0)) {
						// Of these candidates that survived, find the closest one:
						
						distsqr=segment_segment_distsqr(elems[i].pos, elems[i+1].pos,elems[j].pos, elems[j+1].pos, s1p, s2p, s1t, s2t);
						// If the square distance is smaller than the square of thickness*2, they are contact candidates
					
						if (distsqr<tn2sqr) {
						    //Check if two close segments along the wire have their closest point along the wire, ie not really a contact
						  if (i+2<j || (s1t<1. || s2t>0.)) 
						    {	
						      /*    if (sweeps > 900000 && i>20)
							{
							  std::cout<<"Sweep "<<sweeps<<" elem "<<i<<" has contact with "<<j<<". This is contact number "<<elems[i].contacts.size()
								   <<", distsqr: "<<distsqr<<", s1t: "<<s1t<<", s2t: "<<s2t
								   <<std::endl;
							}
						      */
						      //Check if we already have this contact from the last LC rebuild:
						    if (elems[i].ct_idx_map[j]!=1) { 
						    	newcontact.idx=j;
								newcontact.valid=true;
						    	newcontact.distsqr=distsqr;
						    	newcontact.s1t=s1t;
						    	newcontact.s2t=s2t;
						    	newcontact.s1p=s1p;
						    	newcontact.s2p=s2p;
								newcontact.mu=wire_friction_const;
						    	elems[i].contacts.push_back(newcontact);
								elems[i].ct_idx_map[j]=1;
							
						    }

						    //Put the ith element in the global contact list if not already done
						    if (!elems[i].in_ct_list) {
								contact_list.push_back(i);
								elems[i].in_ct_list=true;
							}				
						    
						    }
						}
					    }
					}
					j=linkedcell[j];
				    }
				   
				}
			    }
			i=linkedcell[i];
			}	   
	    	} 
	    }
    } 
    //   std::cout<<"Linked cell completely rebuilt\n";    
    //Now we can search the contact_list:
    double repulsionforcel=0.;
    double thicknesssqr=thickness*thickness;
    
    for (std::vector<int>::iterator it=contact_list.begin(); it<contact_list.end(); ++it) {
	i=*it;
	
    	for (int k=0; k<elems[i].contacts.size(); ++k) {
	    j=elems[i].contacts[k].idx; //j is the index of the contact node
	    if (j<0) { // Just to make sure we're not out of range later
			std::cout<<"Something is wrong with the contact data: Contact partner "<<k<<" of segment "<<i<<" has negative index\n";
			all_done=1.;
			exception_occured=true;
			return;
		}
            if (just_rebuilt)
                distsqr=elems[i].contacts[k].distsqr;
            else
                distsqr=segment_segment_distsqr(elems[i].pos, elems[i+1].pos,elems[j].pos,elems[j+1].pos,s1p,s2p, s1t,s2t);
		if (distsqr<tn2sqr) { //A contact is still here, so update its properties
		    elems[i].contacts[k].distsqr=distsqr;
		    elems[i].contacts[k].valid=true;
		    elems[i].contacts[k].s1t=s1t;
		    elems[i].contacts[k].s2t=s2t;
		    elems[i].contacts[k].s1p=s1p;
		    elems[i].contacts[k].s2p=s2p;
		   
			    calc_normforce_and_friction(i,k);	  
		} else { //The contact is gone out of reach, therefore flag it such and reset the contact spring elongation
		    elems[i].contacts[k].valid=false;
		    elems[i].contacts[k].distsqr=distsqr;
		    elems[i].contacts[k].ctqv=Tvect(0.);
		}
	   // }
	}
    }
    
}

void WireSim::calc_normforce_and_friction(int i, int k) {
    Tvect normalforce;
    wireiter it=elems.begin();
    double normalforcel=0.;
    double vreldvec;
    double fresdvec;
    double frictionl;
    double dvecl;
    double tcom;
    double comdist;
    double maxf;
    double stopforcel;
    int j=elems[i].contacts[k].idx;
    double dist=sqrt(elems[i].contacts[k].distsqr); 
    
    Tvect dvec;
    Tvect ftang,fperp;
    double cpt;
    Tvect vrel;
    Tvect friction(0.);

    Tvect pij=(elems[i].contacts[k].s1p-elems[i].contacts[k].s2p)/dist; //Unit vector from intersect point to intersect point
    //These are the normalized segment vectors
    Tvect s1=(elems[i+1].pos - elems[i].pos)/elems.get_edge(i).length;
    Tvect s2=(elems[j+1].pos - elems[j].pos)/elems.get_edge(j).length;
	//Calculate velocity of contact points (linear ansatz)
	Tvect vpt1 = elems[i].vel - (elems[i].vel - elems[i+1].vel)*elems[i].contacts[k].s1t;
	Tvect vpt2 = elems[j].vel - (elems[j].vel - elems[j+1].vel)*elems[i].contacts[k].s2t;	
	vrel=vpt1-vpt2; //Not yet in tangential plane!!!!    


    if (dist<thickness) { 
		normalforcel=Estar*(thickness-dist);
//		normalforcel=blitz::pow2((dist - thickness))*sphere_force_const;
		
		ContactEnergy +=Estar*blitz::pow2(thickness-dist);
		elems[i].E_ct +=Estar*blitz::pow2(thickness-dist);
		elems[i].E_ct_counter++;
		elems[j].E_ct +=Estar*blitz::pow2(thickness-dist);
		elems[j].E_ct_counter++;

		//sc_poten+=1./3.*blitz::pow3(fabs(dist - thickness))*sphere_force_const;
		CoordinationNumber++;
			  

	} else {
		//If we are outside the friction "smoothness" zone, we should reset the static friction elongation
		if (dist>thickness*1.05) 
			elems[i].contacts[k].ctqv=Tvect(0.);
		//if not (ie we are in the smoothness zone), we just do nothing (but do not reset 
		//the friction spring, just leave it without update).
		//This means that the same static friction will act when we get in contact again
		return;
    }
	normalforce=pij*normalforcel;
//if (sweeps%1000==1)
//	std::cout<<"elem "<<i<<", contact "<<k<<" (elem "<<j<<"): contact normalforce: "<<normalforce;
	//Add some damping to it:
	double dampl=spring_damping_const*blitz::dot(pij,vrel);
	if (normalforcel>dampl) //Only subtract damping if it is smaller than the normalforce
	    normalforce-=pij*dampl;
   
	// Now the friction part: //

	{
	//Now, determine the plane perpendicular to the contact line. We do this in a way to remove numerical problems
	Tvect v1=pij;
	Tvect ort[3];
	ort[0]=Tvect(0.);
	ort[0](0)=1.;
	ort[1]=Tvect(0.);
	ort[1](1)=1.;
	ort[2]=Tvect(0.);
	ort[2](2)=1.;
		
	double val=1000.;
	int ind=1000;
	for(int kk=0; kk< 3; kk++)
	{
	    if(fabs(v1[kk]) < val)
	    {
        	val = fabs(v1[kk]);
	        ind = kk;
	    }
	}
	//Together with v1, we have now a local coord. system, where v1 points in the contact direction (normal)
	//v1,v2 and v3 are all orthonormal vectors and form a right-handed 3-bein
	Tvect v3=blitz::cross(v1,ort[ind]);
	Tvect v2=blitz::cross(v3,v1);
    //Resolve cvvtang along v2 and v3
	double vrelv2l = blitz::dot(vrel,v2);
	double vrelv3l = blitz::dot(vrel,v3);
	Tvect cvvtang=vrelv2l*v2+vrelv3l*v3;
	double cvvtangl2=blitz::dot(cvvtang,cvvtang);

//	double cvvtangl=sqrt(blitz::dot(cvvtang,cvvtang));
//	friction = -cvvtang/cvvtangl*wire_friction_const*normalforcel;
	bool sticky;
	Tvect fdamp(0.);
		
	//if (cvvtangl2<STICKY_MIN_VEL*STICKY_MIN_VEL) {
	    elems[i].contacts[k].ctqv += cvvtang*dt;
	    friction = -STICKY_FORCE_CONST*elems[i].contacts[k].ctqv;
	    //Project it on the perp. plane by removing the normal component:
	    friction -= blitz::dot(friction,v1)*v1;
		
		fdamp=-spring_damping_const*cvvtang;		
		sticky=true;
	//}
		frictionl = blitz::dot(friction,friction);
		if (blitz::dot(fdamp,fdamp)<frictionl)
			friction+=fdamp;
		
	if (frictionl>=blitz::pow2(wire_static_friction_const*normalforcel) || cvvtangl2>STICKY_MIN_VEL*STICKY_MIN_VEL) {
		frictionl=sqrt(frictionl);
		if (elems[i].contacts[k].mu>wire_friction_const)
			elems[i].contacts[k].mu+= (wire_friction_const-wire_static_friction_const)/10.;
		else elems[i].contacts[k].mu=wire_friction_const;
		friction=friction*((elems[i].contacts[k].mu*normalforcel)/frictionl);
		//	friction = it->vel*((it->wmu*normalforcel)/it->vel.mag());
		//	friction -= friction.dot(dvec)*dvec;
		
		fdamp=-spring_damping_const*cvvtang;
		
		frictionl = blitz::dot(friction,friction);
		if (blitz::dot(fdamp,fdamp)<frictionl)
			friction+=fdamp;
		
		elems[i].contacts[k].ctqv=-friction/STICKY_FORCE_CONST;
				
		/*
		 double cvvtangl = sqrt(cvvtangl2);
		if (cvvtangl>1.E-20)
		{
	    friction = -cvvtang/cvvtangl*wire_friction_const*normalforcel;
			friction=-wall_friction_const*normalforcel*vperp/vperplen;
			double frictionl2=blitz::dot(friction,friction);
			//The part of the previous timestep's force along the current friction force direction
			Tvect fprevl=sphere_mass*(blitz::dot(it->ap,vperp)*vperp/(vperpl2)+blitz::dot(it->vel,vperp)*vperp/(vperpl2)/(30.*dt));
			double fprevl2=blitz::dot(fprevl,fprevl);
			//The new friction shouldn't be greater than the previous force
			if (frictionl2>=fprevl2)
			{ //Maximum possible friction force is reached -> rescale friction
				friction*=sqrt(fprevl2/frictionl2);
			    //	std::cout<<"Rescaling friction of element "<<it->ID<<" to "<<friction<<std::endl;
			}
			
		elems[i].contacts[k].ctqv=-friction/STICKY_FORCE_CONST;
			
		}*/
		
		sticky=false;
	}
		if (sticky)
			StaticFrictionEnergy += 0.5*STICKY_FORCE_CONST*blitz::dot(elems[i].contacts[k].ctqv,elems[i].contacts[k].ctqv);

	//Check if both friction coefficients are zero. Only need to check the first, actually:
	if (wire_static_friction_const==0.)
		friction=-tangential_contact_damping_const*cvvtang;
		

    }
	 

	//Now we start applying the normalforce to the segments. Please note that the friction force is always acting only on the
	//center of mass and is not causing any rotation. This just makes it more stable...
	
    cpt=elems[i].contacts[k].s1t;
    //Equally split the translational movement according to the contact point on segment 1:
//Debuggin	
//   elems[i].force+=0.5*(normalforce+friction);
	//elems[i+1].force+=0.5*(normalforce+friction);
elems[i].force+=(1.-cpt)*(normalforce+friction);
elems[i+1].force+=cpt*(normalforce+friction);
    // Okay, now deal with the angular momentum if the force acts not exactly on the center of mass:
 /*   dvec=-elems[i].pos+elems[i+1].pos;
    it=elems.begin()+i;
    dvecl=it.right_edge().length;
    dvec/=dvecl;  //Normlize
    ftang=blitz::dot(dvec, normalforce); //Tangential part of force
    fperp=normalforce-ftang*dvecl; //Perpendicular part of normalforce
    tcom=0.5; //Center of mass is at t=0.5
    comdist=tcom-cpt;
*/
  //Debugging	
//   elems[i].force+=fperp*(comdist)/0.5;
	//elems[i+1].force-=fperp*(comdist)/0.5;

    //Now the second (j) line segment:
    cpt=elems[i].contacts[k].s2t;
    normalforce*=-1.; //Now, normalforce points into the right direction, seen from segment 2
    friction*=-1.;
    //Equally split the translational movement according to the contact point on segment 2:
//Debugging	
//   elems[j].force+=0.5*(normalforce+friction);  // cpt*
	//  elems[j+1].force+=0.5*(normalforce+friction);  // (1-cpt)*

	elems[j].force+=(1.-cpt)*(normalforce+friction);  // cpt*
	elems[j+1].force+=cpt*(normalforce+friction);  // (1-cpt)*
    // Okay, now deal with the angular momentum if the force acts not exactly on the center of momentum:
  /*  dvec=-elems[j].pos+elems[j+1].pos;
    it=elems.begin()+j;
    dvecl=it.right_edge().length; 

    dvec/=dvecl;  //Normlize
    ftang=blitz::dot(dvec, normalforce); //Tangential part of force
    fperp=normalforce-ftang*dvecl; //Perpendicular part of normalforce
    tcom=0.5; //Center of mass is at t=0.5
    comdist=tcom-cpt;
*/
   //Debugging
//    elems[j].force+=fperp*(comdist)/0.5;
//    elems[j+1].force-=fperp*(comdist)/0.5;
 
}


// This is currently done only for a spherical boundary
void WireSim::boundary_force() {
    Tvect ndir;
    Tvect dvec;
    double fresdvec;
    double dvecl;
    double vreldvec;
    double stopforcel;
    double frictionl;
    double normalforcel;
    Tvect normalforce;
    double inner_circle_rad2=blitz::pow2(radius-sphere_radius);
    double ring_xrad=sqrt(blitz::pow2(radius+sphere_radius)-thickness*thickness);
    //Although the node is at distance thickness from the middle line, the contact inside the sphere is a bit closer to the midline due to the sphere
    //curvature:
    double touchpointdist = 2.*sphere_radius*radius/(sphere_radius+radius);
    bool sticky;
    for (wireiter it=elems.begin(); it<elems.end(); ++it) {
      dvec=it->pos-container_center;
      dvecl=blitz::dot(dvec,dvec);	
      if (dvecl>=inner_circle_rad2) {
	Tvect ddvec=-dvec;
	ddvec(0)=0.; //Now, ddvec is in the y-z plane
	//    if (debug) std::cout<<"Node "<<it->ID<<" is contacting the wall (probably). ddvec is "<<ddvec<<std::endl;
	

	double ddist = blitz::dot(ddvec,ddvec);
	if (it->ID!=0 && ddist<=touchpointdist*touchpointdist && it->pos[0]<container_center[0]) { //Test if we are in the "tube" on the left (only insertion channel)
	  if (fabs(dvec[0])<=ring_xrad) {//We are at the entrance door 
	    
	    if (ddist == 0.) {
	      ddvec=Tvect(0.);
	      ddvec(2)=1.;
	      ddist=1.0;
	    } else ddist=sqrt(ddist);
	    
	    //    std::cout<<"This node ("<<it->ID<<") is at the entrance door. pos is"<<it->pos<<", ddist is "<<ddist<<std::endl;
	    Tvect vring=-ddvec/ddist*thickness;
	    
	    //It is always the left channel
	    vring(0)=container_center[0]-ring_xrad;
	    vring(1)+=container_center[1];
	    vring(2)+=container_center[2];
	    
	    //vring is now the closest point on the channel ring of the left or right insertion channel
	    //Now, get the vector from the node to this point on the ring:
	    //  std::cout<<"      vring: "<<vring;
	    vring-=it->pos;
	    ddist=sqrt(blitz::dot(vring,vring));
	    //std::cout<<"      dist(vring-pos): "<<ddist<<std::endl;
	    // Here, we check if the distance from the channel ring is smaller than the sphere radius.
	    // If so, it means that the sphere touches the ring. IMPORTANT: ddist here is not the same as
	    // the ddist checked above.
	    if (ddist<=thickness) {
	      vring/=-ddist;
	      normalforcel=(thickness-ddist)*wall_force_const;
	      normalforce=vring*normalforcel;
	      
	      //Add some damping to it:
	      double dampl=spring_damping_const*blitz::dot(it->vel,vring);
	      if (normalforcel>dampl)  //Only subtract damping if it is smaller than the normalforce!
		normalforce-=dampl*vring;
	      
	      it->force+=normalforce;
	    }
	  } else { //It is in the insertion channel 
	    
	    if (ddist > 0.) {
	      //in direction of ddvec, length: square of distance to channel midline (ddist is here 
	      // still the distance squared!!!
	      ddist=sqrt(ddist);
	      ddvec/=ddist;
	      normalforcel=wall_force_const*ddist;
	      normalforce=normalforcel*ddvec;  // F= k*|d| * d/|d| = k*d 
	      double dampl=spring_damping_const*blitz::dot(it->vel,ddvec);
	      if (normalforcel>dampl)  //Only subtract damping if it is smaller than the normalforce!
		normalforce-=dampl*ddvec;
	      
	      it->force+=normalforce;		      
	    }
	  }
	  if (debug && blitz::dot(normalforce, normalforce)>5)
	    std::cout<<"Element "<<it->ID<<" Channel Boundary Normalforce: "<<normalforce<<std::endl;
	} else {
	  //	if (!it->wctflag) 
	  //    it->wctflag=true;
	  
	  dvecl=sqrt(dvecl);
	  dvec/=dvecl;
	  
	  normalforcel=(dvecl-(radius-sphere_radius))*Estar;
	  normalforce=-dvec*normalforcel;
	  
	  //Add damping to it if it is smaller than the normalforce:
	  double dampl=spring_damping_const*blitz::dot(it->vel,-dvec);
	  if (normalforcel>dampl)  //Only subtract damping if it is smaller than the normalforce!
	    normalforce-=dampl*(-dvec);
	  
	  it->force+=normalforce;
	  
	  Pressure+=normalforcel;
	    
	  Tvect friction(0.);
	  Tvect vperp = it->vel - blitz::dot(dvec,it->vel)*dvec;
	  
	  //Tvect vnorm=blitz::cross(dvec,it->vel); //Normal vector on the plane spanned by position unit vector and velocity
	  //Tvect vperp=blitz::cross(dvec,vnorm); //Tangent to the sphere which is also in the plane spanned by position unit vector and velocity
	  
	  /*	//DEBUGGING:
		Tvect d1=it->vel*dt; //Distance traveled since last timestep
		Tvect dFF=-STICKY_FORCE_CONST*d1;
		dFF -= dFF.dot(dvec)*dvec;
		friction = dFF;
		
		double FfricNorm = friction.mag();
		const double forceNorm = normalforcel;
		
		if (FfricNorm > forceNorm*wall_static_friction_const)
		{ // tangential force greater than static friction -> dynamic
		friction=friction*((wall_friction_const*forceNorm)/FfricNorm);
		
		}
	  */			
	  
	  double vperpl2=blitz::dot(vperp,vperp);
	  //	double vperplen=sqrt(vperpl2);
	  
	  //	friction=vperp/vperlen; //Pre-init... The friction points in the direction of the perpendicular tangent
	  
	  //	vreldvec=blitz::dot(it->vel,vperp); //vel and vperp have negative scalar product by construction
	//	fresdvec=blitz::dot(it->force,dvec);
		frictionl=0.;
		sticky=false;
		Tvect fdamp(0.);
			
	//	if (vperpl2<STICKY_MIN_VEL*STICKY_MIN_VEL) {
			it->wctq+=dt*it->vel;

		    friction=-STICKY_FORCE_CONST*it->wctq;
			friction-=friction.dot(dvec)*dvec; //Subtract normal part
			fdamp=-spring_damping_const*it->vel;
			fdamp-=fdamp.dot(dvec)*dvec;
		    sticky=true;
			//std::cout<<"ID "<<it->ID<<"Sticky, with friction "<<friction<<", wctq: "<<it->wctq<<", dt*it->vel: "<<dt*it->vel<<std::endl;
	//	}
		// Get the length^2 of the sticky friction force:
		frictionl = blitz::dot(friction,friction);
		if (fdamp.mag_sq()<frictionl)
			friction+=fdamp;
			
		if ( frictionl>=blitz::pow2(wall_static_friction_const*normalforcel) || vperpl2>STICKY_MIN_VEL*STICKY_MIN_VEL) {
			frictionl=friction.mag();
			if (it->wmu>wall_friction_const)
				it->wmu+= (wall_friction_const-wall_static_friction_const)/10.;
			else it->wmu=wall_friction_const;
		    friction=friction*((it->wmu*normalforcel)/frictionl);
		//	friction = it->vel*((it->wmu*normalforcel)/it->vel.mag());
		//	friction -= friction.dot(dvec)*dvec;
			
			fdamp=-spring_damping_const*it->vel;
			fdamp-=fdamp.dot(dvec)*dvec;
			frictionl = blitz::dot(friction,friction);
			if (fdamp.mag_sq()<frictionl)
				friction+=fdamp;
			
			/*double frictionl2=blitz::dot(friction,friction);
			//The part of the previous timestep's force along the current friction force direction
			Tvect fprevl=sphere_mass*(blitz::dot(it->ap,vperp)*vperp/(vperpl2)+blitz::dot(it->vel,vperp)*vperp/(vperpl2)/(30.*dt));
			double fprevl2=blitz::dot(fprevl,fprevl);
			//The new friction shouldn't be greater than the previous force
			if (frictionl2>=fprevl2)
			{ //Maximum possible friction force is reached -> rescale friction
				friction*=sqrt(fprevl2/frictionl2);
			    //	std::cout<<"Rescaling friction of element "<<it->ID<<" to "<<friction<<std::endl;
			}
			//DEBUGGING
		//	friction = -0.1*vperp;
			 */
			
		    it->wctq=-friction/STICKY_FORCE_CONST;
		//	}
		    sticky=false;
		/*	std::cout<<"Elem "<<it->ID<<": Unsticky, friction "<<friction;
			std::cout<<"wctq: "<<it->wctq<<", normalf: "<<normalforce<<", it->force "<<it->force;
			std::cout<<", fcoeff: "<<it->wmu<<std::endl;
		 */
		 //		std::cout<<"    fprevl2: "<<fprevl2<<", dotprod: "<<blitz::dot(it->ap,vperp/vperplen)<<std::endl;
		} 
			
			if (sticky) {
				StaticFrictionEnergy += 0.5*STICKY_FORCE_CONST*blitz::dot(it->wctq,it->wctq);
				it->wmu=wall_static_friction_const;
			}
			
		//Check if both friction coefficients are zero. Only need to check the first, actually:
			if (wall_static_friction_const==0.) {
				friction=-tangential_contact_damping_const*vperp;
				it->wctq=Tvect(0.);
				if (tangential_contact_damping_const==0.)
					friction=Tvect(0.);
			}
		it->force+=friction;
		//   if (std::distance(elems.begin(),it)==6)
		//	std::cout<<"\nfrictionl: "<<frictionl<<"\n";   
if (debug && blitz::dot(normalforce, normalforce)>5)
    std::cout<<"Element "<<it->ID<<" Normal Boundary Normalforce: "<<normalforce<<std::endl;
	    }	
	} else { //if (it->wctflag!=false) { //No contact, but we had previously
	    // it->wctflag=false;
	    it->wctq=Tvect(0.);
		it->wmu=wall_friction_const;
	}
    }
}


void WireSim::update_forces() {    
   
    grow();
    wireiter it=elems.begin();
    //Don't forget the first node:
    it->force=Tvect(0.);
    it->Qa.reset();
    it->a=Tvect(0.);
    it->reset_measurements();
    //Now, take care of all the other nodes:
    wireiter it2=it+1;
    for (it=elems.begin(); it<elems.end(); ++it) {
      
      it->reset_measurements();
      it->force=Tvect(0.);
      it->a=Tvect(0.);
      it->Qa.reset();
      if (it<elems.end()-1) 
	{
	  it2=it+1;
	  Tvect vec=it->pos-it2->pos;
	  it.right_edge().length=sqrt(blitz::dot(vec,vec));
	}
    }
	

    computeBeamForces();
    boundary_force();
    update_collision_forces();
    
    treat_special_nodes(); 
	
    //std::cout<<"start: "<<elems[0].vel<<", "<<elems[0].force<<std::endl;
	if (langevinsim) {
        double msqrt=sqrt(sphere_mass);
        for (int i=0; i<total_number_of_nodes-7; i++)
        {
            /*elems[i].force(0) += msqrt*langevin_noise*normalRandom(0.,1.);
            elems[i].force(1) += msqrt*langevin_noise*normalRandom(0.,1.);
            elems[i].force(2) += msqrt*langevin_noise*normalRandom(0.,1.);
             */
            elems[i].force(0) += 2.*msqrt*langevin_noise*(rand01()-0.5);
            elems[i].force(1) += 2.*msqrt*langevin_noise*(rand01()-0.5);
            elems[i].force(2) += 2.*msqrt*langevin_noise*(rand01()-0.5);
        }
        
    }
    for (it=elems.begin(); it<elems.end(); ++it) {
      it->force-=global_friction*it->vel*sphere_mass;  //Dampign relative to mass!
	
	it->a = it->force/sphere_mass;
	if (it<elems.end()-1) {
	  it->waf = acelT(it->Q,it->Qa);
	  
	  it->waf -= global_rotational_friction*IMatmult(it->wf);
	}
	it->finalize_measurements();
    }
    
    if ( debug) {
	std::cout<<"All forces calculated\n";

	for (wireiter it=elems.begin(); it!=elems.end(); ++it) {
	    std::cout<<"particle ID="<<it->ID<<", pos: "<<it->pos<<std::endl;
	    std::cout<<"        vel["<<it->ID<<"] "<<it->vel<<"\n";
	    std::cout<<"        force["<<it->ID<<"] "<<it->force<<"\n";
	    std::cout<<"        waf["<<it->ID<<"] "<<it->waf<<"\n";
	    std::cout<<"        Q: ";
	    it->Q.display();
	    std::cout<<"        Qv: ";
	    it->Qv.display();
	    std::cout<<"        Qa: ";
	    it->Qa.display();
	    
	}
    }
    
}



void WireSim::predictor(void)
{
    //Velocity verlet formulation
    for(wireiter it=elems.begin(); it<elems.end(); ++it)
    {
        it->pos0 = it->pos;
        it->vel0 = it->vel;
        it->ap   = it->a;
        
        it->vel  +=  0.5*it->a*dt;
        it->pos  +=  it->vel*dt;
        
    }
}

void WireSim::corrector(void)
{
    //Velocity verlet formulation
    for (wireiter it=elems.begin(); it<elems.end(); ++it)
    {
        it->vel  +=  0.5*it->a*dt;
        KineticEnergy += 0.5*sphere_mass*blitz::dot(it->vel,it->vel);
    }
}


void WireSim::predictorQ(void)
{
    for (wireiter it=elems.begin(); it<elems.end(); ++it)
    {
        // Verlet angular quaternion integration :
        it->Qv = it->Qv + it->Qaccel*0.5*dt;
        it->Q = it->Q + it->Qv * dt;
        if(fabs(it->Q.norm() -1.0) > 1.e-7)
        {
            it->Q.normalize();
        }
        //Calculate angular velocity in frozen frame, because it is needed for global rotational damping
        Quaternion U = it->Q.conj()*it->Qv; //wf = 2*q.conj()*qv
        it->wf(0)=2.*U.q[1];
        it->wf(1)=2.*U.q[2];
        it->wf(2)=2.*U.q[3];
    }
    if (clamp_free_end_rotation_r3)
    {
        //	std::cout<<	elems[0].wf<<std::endl;
        /*
         elems[0].Qv.reset(); //Clamp the torsion
         elems[0].waf=Tvect(0.);
         elems[0].Qa.reset();
         */
    }
}

void WireSim::correctorQ(void)
{

  if (clamp_free_end_rotation_r3)
    {
      wireiter it = elems.begin();
      Tvect projt = (elems[1].pos-elems[0].pos)/it.right_edge().length;
      double projcomp = blitz::dot(elems[0].w, projt);
    
      //elems[0].w-=projcomp*projt;
      // elems[0].w(2)=0;
    }
  for (wireiter it=elems.begin(); it<elems.end(); ++it)
    {
        Tvect temp = IInvMatmult(it->waf - (it->wf.cross(IMatmult(it->wf))));
        it->Qaccel = (it->Qv*Quaternion(0,it->wf(0),it->wf(1),it->wf(2)))*0.5 +
        it->Q*Quaternion(0.,temp(0), temp(1), temp(2))*0.5;
        it->Qv = it->Qv + it->Qaccel*0.5*dt;
    }
}



void WireSim::predictor_corrector_step() {
    predictor();
    if (!flexible_polymer)
        predictorQ();
    update_forces(); 
  //  treat_special_nodes();

    epssqrmax=0.;
    prec_ID=-1;

    corrector();
    if (!flexible_polymer)
        correctorQ();
 //   treat_special_nodes();

    double avg;
    bool reset_avg=false;

}

inline Tvect WireSim::IMatmult(Tvect v)
{
  Tvect w;
  w(0) = IMat(0,0)*v(0);
  w(1) = IMat(1,1)*v(1);
  w(2) = IMat(2,2)*v(2);
  return w;
}

inline Tvect WireSim::IInvMatmult(Tvect v)
{
  Tvect w;
  w(0) = IInvMat(0,0)*v(0);
  w(1) = IInvMat(1,1)*v(1);
  w(2) = IInvMat(2,2)*v(2);
  return w;
}

// Grows the wire at both ends (ie. the length, total mass etc. will increase)
void WireSim::grow(bool vtk_insertion) {
    Tvect pos(0.);
    double nlength=0.;
  
    wireiter it=elems.end();
    wireiter it2=elems.end()-1;  //This will be the left neighbor

    if (elems[total_number_of_nodes-4].pos[0]+spring_x0/2.>=container_center[0]-radius) { // The third node has entered the circle
	pos=2.*elems[total_number_of_nodes-1].pos-elems[total_number_of_nodes-2].pos;
	
	nlength=sqrt(blitz::dot(elems[total_number_of_nodes-1].pos-elems[total_number_of_nodes-2].pos,
				elems[total_number_of_nodes-1].pos-elems[total_number_of_nodes-2].pos));

	last_insertion_wait_time<<realsweeptime-last_insertion_time;
	last_insertion_time=realsweeptime;
    
	pos(1)=container_center[1]; //Always do this!!
	pos(2)=container_center[2]; //Always do this!!
	std::cout<<"New element "<<total_number_of_nodes<<" at pos: "<<pos<<", last node: "<<elems[total_number_of_nodes-1].pos <<std::endl;
	rebuild_lc=true;
	wire::node newnode;
	wire::edge dummy,newedge;

	newnode.pos=pos;       
	newnode.ID=total_number_of_nodes;
	newedge.length=nlength;

	elems.insert(elems.end(),newnode,dummy);
	it=elems.end()-2;
	setup_beam(it);
	total_number_of_nodes++;   
    }
}


void WireSim::treat_special_nodes() {

    double vreldvec=0.;
    double dvecl;
    double forcemax;
    double faccel;
    int nn=total_number_of_nodes-1;

    if (during_stiffness_measurement) {       
	elems[nn].force(0)=0.;
	elems[nn].vel(0)=0.;
	elems[nn].a=Tvect(0.);
	return;
    }

    if (insertion_force>0.0) {
	elems[nn].force(1)=0.;
	elems[nn].vel(1)=0.;
	elems[nn].force(2)=0.;
	elems[nn].vel(2)=0.;
	elems[nn].force(0)+=insertion_force;
    } else {
	elems[nn].force=Tvect(0.);
	elems[nn].a=Tvect(0.);
	elems[nn].a1=Tvect(0.); 
	elems[nn].a2=Tvect(0.);
	elems[nn].a3=Tvect(0.);
	elems[nn].ap=Tvect(0.);
	elems[nn].vel=Tvect(0.);
	elems[nn].vel(0)=insertion_speed;
	elems[nn].vel0=Tvect(0.);
	elems[nn].vel0(0)=insertion_speed;
    }  
    if (clamp_nozzle_rotation_r3)
      {
	elems[nn-1].Qv.reset(); //Clamp the torsion
	elems[nn-1].waf=Tvect(0.);
	elems[nn-1].Qa.reset();
      }

    elems[nn].Q = elems[nn-1].Q; //Needed since the last element's Q is not yet set by any force or so
    elems[nn].Qv.reset();

    if (initial_deformation ==0. && elems[0].pos[1]<container_center[1]+1.0 && elems[0].pos[2]<container_center[2]+1.0 && sweeps<100000) {
	//elems[0].a = 0.,0.0,0.;
	elems[0].force = 0.;
	elems[0].vel = 0.;
	elems[0].vel(1)= 0.01;
	elems[0].vel(2)=0.005;
    }
	
}


void WireSim::vtk_update() {
}
