#ifndef WIRE_H
#define WIRE_H 1

#include <libconfig.h++>
#include "container.hpp"
#include <iostream>
#include <fstream>
#include "mili/mili.h"
//#include "boost/format.hpp"
//#include <boost/random/variate_generator.hpp>
//#include <boost/random/normal_distribution.hpp>


#include <ext/hash_map>
#include "utils.h"
#include "quaternion.h"
// #include "vtkWireRenderer.h"


class stiffness_measurement {
 public:
  double initial_measurement;
  double last_measurement;
  double k;
  int n;
  int msweeps;
  double min_kin;
  bool measuring;
  bool during_stiffness_measurement;

  void init() {
    measuring=false;
    step=0;
    pts.resize(n);
    pts[0]=initial_measurement;
    ym = k*last_measurement*last_measurement;
    dy = ym/n;
    double x = initial_measurement;
    for (int i = 1; i<n; i++)
      {
	double dx = (-k*x + sqrt(2.*dy*k + k*k*x*x))/(2.*k);
	pts[i] = x+dx;
	x+=dx;
      }
  };
  bool measure_next(double x) { if (x>=pts[step] ) return true; else return false; };
  void advance() { step++; if (step>n-1) { double dx=pts[n-1]-pts[n-2]; pts.push_back(pts[n-1]+dx); n++; } };
  double next_x() { return pts[step]; }
  int step;
   
  std::vector<double> pts;
  double ym;
  double dy;
};

/*
inline alps::IDump& operator >> (alps::IDump& dump, stiffness_measurement& st) {
    
	dump>>st.msweeps>>st.measuring>>st.during_stiffness_measurement>>st.step;
	
    return dump;
}
inline alps::ODump& operator << (alps::ODump& dump, const stiffness_measurement& st) {
	
    dump<<st.msweeps<<st.measuring<<st.during_stiffness_measurement<<st.step;
    return dump;
}
*/
class average {
public:
  average() { total_value=0.; Nvalues=0; maximum=0.;}
  void add(double value) { total_value += value; Nvalues++; if (value>maximum) maximum=value;}
  double avg() { return total_value/Nvalues; }
  double max() { return maximum; }
  void reset() { total_value=0.; Nvalues=0; maximum=0.;}
private:
  double total_value;
  int Nvalues;
  double maximum;
};


namespace wire {
    class contact {
    public:
	contact() {
	    ctqv=Tvect(0.);
	    valid=true;
	    idx=-1;
	}
	int idx;        //index of the first node of the second segment (S2), ie. the second segment will be idx->idx+1
	double distsqr; //Distance^2 to the line segment from contact_ID->right_edge
	double s1t,s2t; //Parameter in (0,1) for the parametrization of the contact point on the line segments S1 and S2
	Tvect s1p,s2p;  //The contact points on the line segments S1 and S2. Can coincide with the end points of the segment
	Tvect ctqv;
	bool valid;
	double mu;
    };

    class node {
    public:
	node() {
	    vel=Tvect(0.);
	    vel0=Tvect(0.);
	    force=Tvect(0.);
	    a=Tvect(0.);
	    w=Tvect(0.);
	    a1=Tvect(0.);
	    a2=Tvect(0.);
	    a3=Tvect(0.);
	    ap=Tvect(0.);
	    waf=Tvect(0.);
	    wf=Tvect(0.);
	    flen=0.;
	    in_ct_list=false;
	    ct_idx_map.resize(0);
	    contacts.resize(0);
	    wctq=Tvect(0.);
	    wctflag=false;
	    E_bend=0.;
	    E_tens=0.;
	    E_tors=0.;
	    E_ct=0.;
	    Q.set(0,0,0);
	    Qv.set(0,0,0);
	    Qa.set(0,0,0);
        Qaccel.set(0,0,0);
	    Qa1.set(0,0,0);
	    Qa2.set(0,0,0);
	    Qa3.set(0,0,0);
	    Qap.set(0,0,0);
	    Q0.set(0,0,0);
	    Qv0.set(0,0,0); 
	}
	void reset_measurements()
	{
	  E_bend=0;
	  E_tens=0;
	  E_tors=0;
	  E_ct=0;
	  E_ct_counter=0;
	}
	void finalize_measurements()
	{
	 
	  if (E_ct_counter>0)
	    E_ct/=E_ct_counter;
	}
  /*
    void load(mili::bostrem<> dump)
    {
        int ctsize;
        dump>>n.ID>>n.pos>>n.vel>>n.a>>n.flen>>n.waf>>n.w
        >>n.pos0>>n.vel0>>n.ap>>n.a1>>n.a2>>n.a3
        >>n.Q>>n.Qv>>n.Qa>>n.Q0>>n.Qv0>>n.Qa1>>n.Qa2>>n.Qa3>>n.Qap>>n.wf
        >>n.wctq>>n.wctflag>>ctsize>>n.u0;
        
        n.contacts.resize(ctsize);
        //	std::cout<<"Node "<<n.ID<<" dump has "<<ctsize<<" contacts\n";
        if (ctsize>0)
            for (int i=0; i<ctsize; ++i) {
                dump>>n.contacts[i].idx>>n.contacts[i].ctqv>>n.contacts[i].valid>>n.contacts[i].mu;
                //	std::cout<<"contact "<<i<<": idx "<<n.contacts[i].idx<<", valid: "<<n.contacts[i].valid<<"\n";
            }

    
    }
        
        void save(mili::bostrem<> &dump) {
            int ctsize = n.contacts.size();
            dump<<n.ID<<n.pos<<n.vel<<n.a<<n.flen<<n.waf<<n.w
            <<n.pos0<<n.vel0<<n.ap<<n.a1<<n.a2<<n.a3
            <<n.Q<<n.Qv<<n.Qa<<n.Q0<<n.Qv0<<n.Qa1<<n.Qa2<<n.Qa3<<n.Qap<<n.wf
            <<n.wctq<<n.wctflag<<ctsize<<n.u0;
            //	std::cout<<"Node "<<n.ID<<": dumping "<<ctsize<<" contacts\n";
            if (ctsize>0) {
                for (int i=0; i<n.contacts.size(); ++i) {
                    dump<<n.contacts[i].idx<<n.contacts[i].ctqv<<n.contacts[i].valid<<n.contacts[i].mu;
                    //	std::cout<<"contact "<<i<<": idx "<<n.contacts[i].idx<<", valid: "<<n.contacts[i].valid<<"\n";
                }
            }
            
        }
*/
	int ID; 
	Tvect pos;
	Tvect vel;
	Tvect force;
	Tvect a;
	Tvect pos0;
	Tvect vel0;
	Tvect a1,a2,a3,ap; //Used for predictor-corrector
	Tvect waf; //Angular acceleration in frozen system
	double flen;
	double wmu;
	std::vector<contact> contacts;
	__gnu_cxx::hash_map<int, int> ct_idx_map;
	Tvect wctq;
	bool wctflag;
	double E_bend;
	double E_tens;
	double E_tors;
	double E_ct;
	int E_ct_counter;

	bool in_ct_list;
	Quaternion Q, Qv, Qa, Qaccel; //Quaternion to describe the rotation, velocity and acceleration
	Quaternion Qa1,Qa2,Qa3,Qap,Q0,Qv0; //Quaternions needed for predictor-corrector
	Tvect w;
	Tvect wf;
	Tvect u0;
    };
    
    
    class edge {
    public:
	edge() {
	   
	}
	int ID;
	double length;
    };
      
}

class Parameters {
 public:
	Parameters(std::string filename);
	Parameters() {};
	void setup(std::string filename);
	libconfig::Config cfg;
};


class WireSim
{
 public: 
    WireSim() { std::cout<<"ERROR\n"; exit(0); };
    WireSim(std::string filename);
    ~WireSim() {this->halt();};
    void read_params();
    bool dumping;
    void initialize();
    std::vector<Tvect> dump_contact_forces(std::vector<Tvect> positions);
    void dostep();
    std::vector<Tvect> dump_contact_list(std::vector<Tvect> positions);
    std::vector<wire::node> get_all_nodes() {return elems.get_nodes(); }
  
    void reset_nodes();
    void save() const;
    //void save_to_file(const boost::filesystem::path&) const;
    void load();
    double work_done() const;
    void start();
    void halt();
    //bool change_parameter(const std::string& name, const alps::StringValue& value);
    typedef wirecontainer<wire::node,wire::edge>::wireiterator wireiter;
    void vtk_update();
    void setup_beam(wireiter &it);
    void output_vtk(std::string filename);
    Parameters p;
 private:
    std::string outfile;
    std::string infile;

    inline Tvect IMatmult(Tvect v);
    inline Tvect IInvMatmult(Tvect v);

	/**
	 * Measurement variables that are averaged over time
	 */
	average av_BeamBendEnergy;
	average av_BeamCouplEnergy;
	average av_BeamElongEnergy;
	average av_BeamTorsEnergy;
	average av_ContactEnergy;
	average av_InsertionForce;
	average av_Pressure;
	average av_CoordinationNumber;
	average av_StaticFrictionEnergy;
	average av_KineticEnergy;
        
	int rebuild_lc_block;
        int seed;	
	double Pressure;
	double ContactEnergy;
	double BeamBendEnergy;
	double BeamElongEnergy;
	double BeamTorsEnergy;
	double BeamCouplEnergy;
	double StaticFrictionEnergy;
	double KineticEnergy;
	double CoordinationNumber;

	stiffness_measurement st;
	std::string st_measure_file;
	std::string measure_file;

	void reset_all_measurement_averages();
	void reset_all_measurement_variables();

	void computeBeamForces();

	double STICKY_FORCE_CONST;
	average av_epssqrmax;
	double epssqrmax; /** Maximum corrctor norm^2 */
	double prec_minimum_tolsqr;
	double prec_maximum_tolsqr;
	double prec_maximum_dt;
	double prec_minimum_dt;
	int prec_last_increase_sweep; /** Last iteration when we increased the time step */
	int prec_ID;
	bool fixed_dt;
	
	double container_volume;  //Both needed for packing_density calculation
	double wire_volume_density;
	double packing_ratio;

	//   vtkWireRenderer vtk_wire;
	// pthread_t vtk_thread;
	// pthread_attr_t attr;
    double sphere_radius;
    std::vector<int> contact_list;
    std::vector<Tvect> dumping_ct_list;
    double all_done;
    MovAvg<double> last_insertion_wait_time;
    double last_insertion_time;
    double insertion_force_start;
    double insertion_force_end;
    double insertion_force_incr_wait_factor;
    double insertion_force_incr_min_wait_time;
    double insertion_force_incr_max_wait_time;
    double insertion_force_incr_factor;
    bool force_increment_sim;
    bool force_increment_phase;
    bool first_wall_contact;
    bool exception_occured;
    bool insert_asymmetric;
    std::ofstream imagefile;
    std::ofstream ctfile;
    std::ofstream measurefile;
    std::string file_prefix;
    double initial_deformation;
 
    double torsion_plastic_factor;
    double global_friction;
    double global_rotational_friction;
    int plot_block;
    int plot_position_block;
    int measure_block;
    bool debug;
    int MAXID;
    int total_number_of_nodes;
    Tvect container_center;
    double radius; //Radius squared of the boundary circle
    double radius2;
    double cyl_radius;
    double cyl_radius2;
    double cyl_length;
    double cyl_length2;
    double thickness;
    double channel_radius;
    double wall_force_const;
    double res_max_angle;
    double res_min_angle;
    double L_x;   // Simulation box width 
    double L_y;   // Simulation box height
    double L_z;   // Simulation box height
    double high_res_length;
    double low_res_length;
    double max_angle0;
    double vmaxvmax; // Maximum velocity squared (is used for linked cell size calculation
    double deform_angle_thresh; //The angle above which plastic deformation occurs
    int pcorder;     // Predictor-corrector order
    wirecontainer<wire::node,wire::edge> elems;
    
    double emodul;
	double Estar; // Effective modulus for cylinder contacts
    double sphere_inertia; //Inertia for the nodes (all the same)
    double sphere_mass;
    double spring_x0;
    double position_damping_const;
    double bending_force_const;
    double torsion_force_const;
    double spring_force_const;
    double bending_damping_const;
    double torsion_damping_const;
    double spring_damping_const;
    double coupling_const;
    double tangential_contact_damping_const;
    double u0x_strain, u0y_strain, u0z_strain;
    double last_u0z_strain;
    double LM_u0xcurrent;
    int LM_segmentctr;
    int LM_nsegments;
    
    bool langevinsim;
    double temperature;
    double langevin_gamma;
    double langevin_noise;

    double LM_avgloopsize;
    bool LM_loopmodel;
    double u0x_sigma;
    double u0y_sigma;
    double u0z_sigma;
    bool eqmodel;
    bool u0z_incr_random;
    bool hamiltonwalk;
    double hamiltonu0;
    int hamilton_ctr;
    double qyieldthresh;
    double strain_hardening_factor;
    double wall_friction_const;
    double wire_friction_const;
    double wire_static_friction_const;
    double wall_static_friction_const;
 
    bool block_deform; //Used to block plastic deformation while relaxing resolution changes
    int block_updates; 
    int block_sweeps;
    Tmat IMat, IInvMat;
    Matrix44 B1,B2,B3, B01,B02,B03;

    double dt;
    int sweeps;
    int total_sweeps;
    void grow(bool vtk_insertion=true);  
    void treat_special_nodes();
    void setup_matrices();
  
    void predictor_corrector_step();
    void predictor();
    void predictorQ();
    void corrector();
    void correctorQ();
    
    void do_measurements(std::string filename, bool write_measurements);
    
    void update_forces();
    void boundary_force();
    double segment_segment_distsqr(Tvect S1P0, Tvect S1P1, Tvect S2P0, Tvect S2P1, Tvect &s1p, Tvect &s2p, double &s1t, double &s2t);
    
    bool strain_hardening(double &f, double &angle, double &angleP, double &amin, double &amax);
    
    void update_collision_forces();
    void calc_normforce_and_friction(int i, int k);
    void draw_nodes(); 
    std::vector<int> head;
    std::vector<int> linkedcell;
    double rc;
    double rcrc;
    int L_cx;
    double largest_poten;
    int L_cy;
    int L_cz;
    bool rebuild_lc;
    double insertion_speed;
    std::deque<double> MovAvForceLeft; //Last external forces at left side. Used for moving average calculation
    std::deque<double> MovAvForceRight; //Last external forces at right side. Used for moving average calculation
    int NMovAv; // Moving average length
    double MovAvForceLeftAvg;
    double MovAvForceRightAvg;
    int measurements_done;
    double realsweeptime;
    bool measure_crit_len;
    double tn2;
    double tn2sqr; 
    bool restarted;
    bool use_channel_friction;
    double insertion_force;
    double tnextcol;
    double lambda;
    bool soft_core_potential;
    bool dense_packing;
    double insertion_speed_factor;
    double max_packing;
    int old_total_number_of_nodes;
    double pressure_poten;
    double sphere_force_const;
    bool xwin_output;
    int max_length;
    bool measure_stiffness;
    bool during_stiffness_measurement;
    double stiffness_left;
    double stiffness_right;
    int measure_stiffness_block;
    int measure_stiffness_relaxation_block;
    bool clamp_nozzle_rotation_r3;
    bool clamp_free_end_rotation_r3;
    bool flexible_polymer;
};



inline mili::bistream<>& operator >> (mili::bistream<>& dump, Tvect& x) {
    for (int i=0;i<3;++i)
	dump>>x(i);
    return dump;
}

inline mili::bostream<>& operator << (mili::bostream<>& dump, const Tvect& x) {
    for (int i=0;i<3;++i)
	dump<<x[i];
    return dump;
}

inline mili::bistream<>& operator >> (mili::bistream<>& dump, Quaternion& q) {
    double q1,q2,q3,q4;
    dump>>q1>>q2>>q3>>q4;
    q.set(q1,q2,q3,q4);
    return dump;
}
inline mili::bostream<>& operator << (mili::bostream<>& dump, const Quaternion& q) {
    dump<<q.q[0]<<q.q[1]<<q.q[2]<<q.q[3];
    return dump;
}

inline mili::bistream<>& operator >> (mili::bistream<>& dump, wire::node& n) {
	int ctsize;
    dump>>n.ID>>n.pos>>n.vel>>n.a>>n.flen>>n.waf>>n.w
	>>n.pos0>>n.vel0>>n.ap>>n.a1>>n.a2>>n.a3
	>>n.Q>>n.Qv>>n.Qa>>n.Q0>>n.Qv0>>n.Qaccel>>n.Qa2>>n.Qa3>>n.Qap>>n.wf
	>>n.wctq>>n.wctflag>>ctsize>>n.u0;

    n.contacts.resize(ctsize);
//	std::cout<<"Node "<<n.ID<<" dump has "<<ctsize<<" contacts\n";
    if (ctsize>0)
	for (int i=0; i<ctsize; ++i) {
	    dump>>n.contacts[i].idx>>n.contacts[i].ctqv>>n.contacts[i].valid>>n.contacts[i].mu;
	    //	std::cout<<"contact "<<i<<": idx "<<n.contacts[i].idx<<", valid: "<<n.contacts[i].valid<<"\n";
	}
	
    return dump;
}
inline mili::bostream<>& operator << (mili::bostream<>& dump, const wire::node& n) {
    int ctsize = n.contacts.size();
    dump<<n.ID<<n.pos<<n.vel<<n.a<<n.flen<<n.waf<<n.w
	<<n.pos0<<n.vel0<<n.ap<<n.a1<<n.a2<<n.a3
	<<n.Q<<n.Qv<<n.Qa<<n.Q0<<n.Qv0<<n.Qaccel<<n.Qa2<<n.Qa3<<n.Qap<<n.wf
	<<n.wctq<<n.wctflag<<ctsize<<n.u0;
//	std::cout<<"Node "<<n.ID<<": dumping "<<ctsize<<" contacts\n";
    if (ctsize>0) {
	for (int i=0; i<n.contacts.size(); ++i) {
	    dump<<n.contacts[i].idx<<n.contacts[i].ctqv<<n.contacts[i].valid<<n.contacts[i].mu;
	    //	std::cout<<"contact "<<i<<": idx "<<n.contacts[i].idx<<", valid: "<<n.contacts[i].valid<<"\n";
	}
	}
	return dump;
}

inline mili::bistream<>& operator >> (mili::bistream<>& dump, wire::edge& e) {
    
    dump>>e.ID>>e.length;
    return dump;
}
inline mili::bostream<>& operator << (mili::bostream<>& dump, const wire::edge& e) {
    
  dump<<e.ID<<e.length;
    return dump;
}

#endif

